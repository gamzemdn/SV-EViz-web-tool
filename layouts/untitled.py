 
    # --- List output files
    if not out_dir or not os.path.isdir(out_dir):
        return status_div, html.Ul([html.Li("Output directory not found.")]), dash.no_update

    try:
        names = sorted(os.listdir(out_dir))
        if summary_data:
            df_summary = pd.DataFrame(list(summary_data.items()), columns=["Metric", "Value"])
            metrics_figs = generate_truvari_visuals(df_summary)
        else:
            metrics_figs = html.Div("❌ No summary.json data found.", style={"color": "red"})
      
        return status_div, html.Ul([]), metrics_figs
    except Exception as e:
        return status_div, html.Ul([html.Li(f"Listing error: {e}")]), dash.no_update
















        @app.callback(
    Output('truvari-visual-output', 'children'),
    Input('truvari-uploaded-vcf-paths', 'data')
)
def update_truvari_visuals(vcf_paths):
    if not vcf_paths:
        return dash.no_update

    # Dosya yolları
    tp_comp_path = vcf_paths.get("tp_comp")
    tp_base_path = vcf_paths.get("tp_base")
    fp_path = vcf_paths.get("fp")
    fn_path = vcf_paths.get("fn")

    plots = []

    # TP-comp
    df_tp = load_vcf_dataframe(tp_comp_path)
    plots.append(html.Hr())
    plots.append(html.H4("✅ TP-comp Manhattan"))
    plots.append(plot_manhattan(df_tp))

    # FP
    df_fp = load_vcf_dataframe(fp_path)
    plots.append(html.Hr())
    plots.append(html.H4("❌ FP Clustergram"))
    plots.append(plot_clustergram(df_fp, selected_chroms=["chr1", "chr2", "chr3"]))

    # FN
    df_fn = load_vcf_dataframe(fn_path)
    plots.append(html.Hr())
    plots.append(html.H4("⚠️ FN Sankey"))
    plots.append(plot_sankey(df_fn))

    # TP-base örneği: Circos
    circos_json_path = "tp_base_circos.json"
    vcf_to_circos_json(tp_base_path, circos_json_path, bin_size=1000000)
    plots.append(html.Hr())
    plots.append(html.H4("🧬 TP-base Circos"))
    plots.append(plot_circos(circos_json_path, "DEL"))

    return html.Div(plots)



















    
    def _section(label, vcf_path):
        if not vcf_path or not os.path.exists(vcf_path):
            return None
        try:
            df = load_vcf_dataframe(vcf_path, "caller")  # truvari çıktıları standart VCF
            if df is None or df.empty:
                return html.Div(f"⚠️ {label}: no data.")
            # Basit görseller (bar, radar, histogram, violin vs.)
            basic_figs = plot_vcf_data(df)

            # Sankey
            sankey_html = plot_sankey(df)

            # Clustergram (mevcut kromlar otomatik)
            valid_chroms = [f"chr{i}" for i in range(1,23)] + ["chrX","chrY"]
            df_chr = df.copy()
            df_chr["CHROM"] = df_chr["CHROM"].astype(str).apply(lambda x: f"chr{x}" if not str(x).startswith("chr") else x)
            available_chroms = [c for c in valid_chroms if c in df_chr["CHROM"].unique()]
            cluster_html = plot_clustergram(df_chr, available_chroms) if available_chroms else html.Div("⚠️ No chromosomes for clustergram.")

            # Circos JSON + Circos plot
            json_path = vcf_path + "_circos.json"

            plain_vcf_for_circos = _ensure_plain_vcf_for_circos(vcf_path)
            vcf_to_circos_json(plain_vcf_for_circos, json_path, "truvari")
            
            with open(json_path, "r") as f:
                circos_data = json.load(f)
            
            # SV type renk haritası çıkar
            svtype_colors = {track["name"]: track.get("color", "#95a5a6") for track in circos_data.get("tracks", [])}
            
            # SVTYPE listesi (örnek: df içinde varsa alabilirsin)
            svtypes_list = df["SVTYPE"].dropna().unique().tolist() if "SVTYPE" in df.columns else []
            
            circos_html = plot_circos("histogram", svtypes_list, circos_data, svtype_colors)
            # Manhattan (QUAL tabanlı - fonksiyon içinde hesaplanıyor)
            try:
                manhattan_html = plot_manhattan(df, svtypes_list, 6)
            except Exception as me:
                manhattan_html = html.Div(f"❌ Manhattan error: {me}")

            return html.Div([
                html.H4(label, style={'marginTop': '15px'}),
                html.Div(basic_figs),
                html.Br(),
                sankey_html,
                html.Br(),
                html.H4("Clustergram"),
                cluster_html,
                html.Br(),
                html.H4("Circos"),
                circos_html,
                html.Br(),
                html.H4("Manhattan"),
                manhattan_html
            ])
        except Exception as e:
            return html.Div(f"❌ Visualization error for {label}: {e}", style={"color": "red"})

    visuals_sections = []
    for label, key in [("✅ TP-comp", "tp_comp"),
                       ("✅ TP-base", "tp_base"),
                       ("❌ FP", "fp"),
                       ("⚠️ FN", "fn")]:
        sec = _section(label, tru_paths.get(key))
        if sec:
            visuals_sections += [sec, html.Hr()]

    visuals_block = html.Div(visuals_sections) if visuals_sections else html.Div("⚠️ No visualization data found.")

    # Dosya listesi istenmiyor -> boş
    return status_div, html.Ul([]), metrics_preview, visuals_block




























    sankey_df = data.dropna(subset=["CHROM", "SVTYPE"]).copy()
    sankey_df["CHROM"] = "chr" + sankey_df["CHROM"].astype(str)
    sankey_df["CHROM_SRC"] = sankey_df["CHROM"] + " (" + sankey_df["__SOURCE__"] + ")"
    sankey_counts = sankey_df.groupby(["CHROM_SRC", "SVTYPE"]).size().reset_index(name="Count")

    chrom_labels = sankey_counts["CHROM_SRC"].unique().tolist()
    svtype_labels = sankey_counts["SVTYPE"].unique().tolist()
    all_labels = chrom_labels + svtype_labels
    label_map = {label: i for i, label in enumerate(all_labels)}

    sankey_fig = go.Figure(go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=all_labels,
            color=["#aaaaaa"] * len(chrom_labels) + [px.colors.qualitative.Vivid[i % 10] for i in range(len(svtype_labels))]
        ),
        link=dict(
            source=sankey_counts["CHROM_SRC"].map(label_map),
            target=sankey_counts["SVTYPE"].map(label_map),
            value=sankey_counts["Count"],
            color="rgba(160,160,160,0.4)"
        )
    ))
    sankey_fig.update_layout(
        title=f"{title} — Chromosome ➜ SVTYPE Sankey",
        font_size=12,
        height=600,
        margin=dict(l=20, r=20, t=60, b=30)
    )





















        clust_df = data.copy()
        pivot = clust_df.pivot_table(index="CHROM", columns="SVTYPE", aggfunc="size", fill_value=0)
        
        # 1) Figure'ı üret
        clustergram_fig = Clustergram(
            data=pivot.values,
            row_labels=pivot.index.tolist(),      # örn: ["chr1", ...]
            column_labels=pivot.columns.tolist(), # örn: ["DEL","INS",...]
            color_map="RdBu",
            height=700,
            width=900,
            # etiketlere yer açmak için oranları büyüt (heatmap:dendrogram)
            display_ratio=[0.85, 0.15],
            # sadece colorbar etiketini gizlemek istersen "color" bırak;
            # tamamen hiç gizleme istemiyorsan [] kullan
            hidden_labels=[],
            color_threshold={"row": 250, "col": 700}
        )
        
        # 2) Tüm eksenlerde etiketleri açık hale getir
        for k, v in clustergram_fig.layout.items():
            # layout'taki xaxis*, yaxis* alanlarını döngüyle yakala
            if isinstance(v, dict):
                continue
        for ax_name in list(clustergram_fig.layout):
            if ax_name.startswith("xaxis") or ax_name.startswith("yaxis"):
                clustergram_fig.layout[ax_name].update(showticklabels=True)
        
        # 3) Biraz margin ver (etiketler taşmasın)
        clustergram_fig.update_layout(margin=dict(l=80, r=40, t=40, b=80))
        
        # 4) Render et
        return html.Div([
            html.H4(f"{title} — Advanced Graphs", style={'marginTop': '15px'}),
            dcc.Graph(figure=sankey_fig),
            dcc.Graph(figure=clustergram_fig),
        ])


















       tp_base_vcf = next((p for p, lab in pairs if lab == "TP-BASE"), None)
        tp_comp_vcf = next((p for p, lab in pairs if lab == "TP-COMP"), None)
        
        circos_div = html.Div("⚠️ Circos: TP-BASE / TP-COMP not found.")
        if tp_base_vcf and tp_comp_vcf and os.path.exists(tp_base_vcf) and os.path.exists(tp_comp_vcf):
            import tempfile, uuid
            tmp_json = os.path.join(tempfile.gettempdir(), f"circos_truvari_{uuid.uuid4().hex}.json")
        
            # JSON üret
            vcf_to_circos_json_truvari(
                tp_base_vcf, tp_comp_vcf,
                tmp_json,
                svtype_color_map=svtype_color_map,
                chrom_color_map=chrom_color_map,
                hg38_lengths=hg38_lengths,
                hg19_lengths=hg19_lengths,
                bin_size=10_000_000
            )
        
            with open(tmp_json) as f:
                circos_data = json.load(f)
        
            layout = circos_data.get("layout", [])
            layout_ids = {str(b.get("id")) for b in layout}
            all_tracks = circos_data.get("tracks", [])
        
            # Her kombinasyona (örn. "DEL (TP-BASE)") farklı renk
            big_palette = (
                px.colors.qualitative.Vivid
                + px.colors.qualitative.Set2
                + px.colors.qualitative.Plotly
                + px.colors.qualitative.Safe
                + px.colors.qualitative.Bold
                + px.colors.qualitative.Pastel
            )
            def combo_color(name: str) -> str:
                return big_palette[abs(hash(name)) % len(big_palette)]
        
            # BASE ve COMP için ayrı bin listeleri
            base_bins, comp_bins = [], []
        
            for t in all_tracks:
                if str(t.get("type", "")).lower() != "histogram":
                    continue
                name = t.get("name", "")
                lname = name.lower()
                is_base = ("tp-base" in lname) or name.endswith(" (TP-BASE)")
                is_comp = ("tp-comp" in lname) or name.endswith(" (TP-COMP)")
                if not (is_base or is_comp):
                    continue
        
                col = combo_color(name)
                for item in t.get("data", []):
                    blk = str(item.get("block_id"))
                    if blk not in layout_ids:
                        continue
                    try:
                        s = int(item.get("start", 0))
                        e = int(item.get("end", 0))
                        v = int(item.get("value", 0))
                    except Exception:
                        continue
                    if e <= s or v < 0:
                        continue
                    datum = {"block_id": blk, "start": s, "end": e, "value": v, "color": col}
                    (base_bins if is_base else comp_bins).append(datum)
        
            tracks = []
            # ✅ TP-BASE = içteki bant (ideogramın içinde)
            if base_bins:
                tracks.append({
                    "type": "HISTOGRAM",              # bazı sürümlerde büyük harf şart
                    "name": "TP-BASE (all SVTYPE combos)",
                    "data": base_bins,
                    "innerRadius": 160,               # ideogram iç yarıçapından küçük → içeride
                    "outerRadius": 195,
                    "config": {                       # geriye dönük uyumluluk
                        "innerRadius": 160,
                        "outerRadius": 195,
                        "opacity": 1.0,
                        "strokeWidth": 0
                    }
                })
        
            # ✅ TP-COMP = onun üstünde ikinci iç bant
            if comp_bins:
                tracks.append({
                    "type": "HISTOGRAM",
                    "name": "TP-COMP (all SVTYPE combos)",
                    "data": comp_bins,
                    "innerRadius": 200,
                    "outerRadius": 235,
                    "config": {
                        "innerRadius": 200,
                        "outerRadius": 235,
                        "opacity": 1.0,
                        "strokeWidth": 0
                    }
                })
        
            if not tracks:
                circos_div = html.Div("⚠️ Circos: no histogram data available.")
            else:
                # Kombinasyon legend’i
                combo_names = sorted({
                    t.get("name", "") for t in all_tracks
                    if str(t.get("type", "")).lower() == "histogram"
                })
                legend = html.Div([
                    html.Div([
                        html.Span(style={
                            'backgroundColor': combo_color(lbl),
                            'display': 'inline-block', 'width': '14px', 'height': '14px',
                            'marginRight': '6px', 'borderRadius': '3px'
                        }),
                        html.Span(lbl, style={'marginRight': '12px'})
                    ], style={'display': 'inline-block', 'marginBottom': '4px'})
                    for lbl in combo_names
                ], style={'padding': '6px 0'})
        
                circos_div = html.Div([
                    html.Hr(),
                    html.H4("Circos — TP-BASE (inner band) & TP-COMP (inner band above)"),
                    legend,
                    dashbio.Circos(
                        id="circos-truvari",
                        layout=layout,
                        tracks=tracks,
                        config={
                            # ideogram dış halkası; histogram bantları bunun İÇİNDE kalacak
                            "innerRadius": 260,        # ideogram iç yarıçapı
                            "outerRadius": 320,        # ideogram dış yarıçapı
                            "ticks": {"display": True, "spacing": 10_000_000, "labelSpacing": 5},
                            "labelLayout": {"spacing": 60, "radialOffset": 90}
                        },
                        size=800
                    )
                ])





































































@app.callback(
    [Output("tru-status", "children"),
     Output("truvari-output-list", "children"),
     Output("truvari-metrics-preview", "children"),
     Output("truvari-visual-output", "children")
    ],  
    Input("tru-run-btn", "n_clicks"),
    State("base_choice", "value"),  # <-- yeni eklendi
    State("tru-base-upload", "filename"),
    State("tru-base-upload", "contents"),
    State("tru-comp-upload", "filename"),
    State("tru-comp-upload", "contents"),
    State("tru-ref-selector", "value"),
    State("tru-param-refdist", "value"),
    State("tru-param-pctsim", "value"),
    State("tru-param-pctsize", "value"),
    State("tru-param-pctovl", "value"),
    State("tru-param-sizemin", "value"),
    State("tru-param-sizefilt", "value"),
    State("tru-param-sizemax", "value"),
    State("tru-flag-typeignore", "value"),
    State("tru-flag-uselev", "value"),
    State("tru-flag-gtcomp", "value"),
    State("tru-flag-passonly", "value"),
    State("tru-flag-multimatch", "value"),
    prevent_initial_call=True
)
def run_truvari_bench(n_clicks,base_choice,
                      base_fname, base_contents,
                      comp_fname, comp_contents,
                      reference_value,
                      refdist, pctsim, pctsize, pctovl, sizemin, sizefilt, sizemax,
                      typeignore, uselev, gtcomp, passonly, multimatch):
    if not n_clicks:
        return dash.no_update, dash.no_update, dash.no_update, dash.no_update

    # --- Validate inputs
    # --- Handle base_vcf based on base_choice
    if base_choice == "univar":
        base_vcf_path = os.path.join(UPLOAD_DIRECTORY, "univar.vcf")
        if not os.path.exists(base_vcf_path):
            return html.Div("❌ Univar SV catalog (univar.vcf) not found."), html.Ul([html.Li("No output.")]), dash.no_update, dash.no_update
    elif base_choice == "custom":
        if not base_contents or not base_fname:
            return html.Div("❌ Please upload the BASE VCF."), html.Ul([html.Li("No output.")]), dash.no_update, dash.no_update
        try:
            base_vcf_path = save_uploaded_file(base_fname, base_contents)
        except Exception as e:
            return html.Div(f"❌ BASE file save error: {e}"), html.Ul([html.Li("No output.")]), dash.no_update, dash.no_update
    else:
        return html.Div("❌ Invalid base choice."), html.Ul([html.Li("No output.")]), dash.no_update, dash.no_update
        
    if not comp_contents or not comp_fname:
        return html.Div("❌ Please upload the COMPARISON VCF."), html.Ul([html.Li("No output.")]), dash.no_update, dash.no_update
    if not reference_value:
        return html.Div("❌ Please select a reference (.fa)."), html.Ul([html.Li("No output.")]), dash.no_update, dash.no_update

    # --- Save uploaded VCFs
    try:
        base_vcf_path = save_uploaded_file(base_fname, base_contents)
        comp_vcf_path = save_uploaded_file(comp_fname, comp_contents)
    except Exception as e:
        return html.Div(f"❌ File save error: {e}"), html.Ul([html.Li("No output.")]), dash.no_update, dash.no_update

    # --- Absolute reference path (layout stores relative 'uploaded_files/..')
    reference_fa = os.path.abspath(reference_value)

    # --- Build Truvari params
    def _maybe(v):  # pass only if user provided
        return v if v not in (None, "") else None

    params = {}
    # flags
    if bool(typeignore): params["--typeignore"] = True
    if bool(uselev):     params["--use-lev"]    = True
    if bool(gtcomp):     params["--gtcomp"]     = True
    if bool(passonly):   params["--passonly"]   = True
    if bool(multimatch): params["--multimatch"] = True
    # keyed values
    if _maybe(refdist)  is not None: params["--refdist"]  = refdist
    if _maybe(pctsim)   is not None: params["--pctsim"]   = pctsim
    if _maybe(pctsize)  is not None: params["--pctsize"]  = pctsize
    if _maybe(pctovl)   is not None: params["--pctovl"]   = pctovl
    if _maybe(sizemin)  is not None: params["--sizemin"]  = sizemin
    if _maybe(sizefilt) is not None: params["--sizefilt"] = sizefilt
    if _maybe(sizemax)  is not None: params["--sizemax"]  = sizemax

    # --- Run pipeline
    status_div, tru_paths = run_truvari_pipeline(
        base_choice,
        base_vcf_uploaded=base_vcf_path,
        comp_vcf=comp_vcf_path,
        reference_fa=reference_fa,
        params=params
    )
    # Pipeline başarısızsa
    if not isinstance(tru_paths, dict):
        return status_div, html.Ul([]), dash.no_update, dash.no_update

    # --- summary.json'u oku ve metrik görsellerini oluştur
    metrics_preview = dash.no_update
    try:
        summary_path = tru_paths.get("summary")
        if summary_path and os.path.exists(summary_path):
            with open(summary_path) as f:
                summary_data = json.load(f)
            # gt_matrix gibi nested yapıları dışla
            df_summary = pd.DataFrame(list(summary_data.items()), columns=["Metric", "Value"])
            df_summary = df_summary[~df_summary["Value"].apply(lambda x: isinstance(x, (dict, list)) or x is None)].copy()
            metrics_preview = generate_truvari_visuals(df_summary)
        else:
            metrics_preview = html.Div("❌ summary.json not found.", style={"color": "red"})
    except Exception as e:
        metrics_preview = html.Div(f"❌ Metrics render error: {e}", style={"color": "red"})

  # --- TRUVARI GÖRSELLERİ (tp-comp, tp-base, fp, fn)
    import plotly.express as px
    from dash import dcc
    
    def _load_and_tag(vcf_path, tag):
        if not vcf_path or not os.path.exists(vcf_path):
            return None
    
        df = load_vcf_dataframe(vcf_path, "caller")  # Truvari çıktı VCF'leri
        if df is None or df.empty:
            return None
    
        # ✅ CHROM normalize
        df["CHROM"] = df["CHROM"].astype(str).str.strip().str.replace("^chr", "", regex=True)
    
        # ✅ SVTYPE'i HER ZAMAN INFO'dan yeniden çıkar (ID/SMAP karışıklığını engeller)
        df["SVTYPE"] = df["INFO"].astype(str).str.extract(r'(?:^|;)SVTYPE=([^;]+)', expand=False)
    
        # (İstersen BND'leri tamamen dışla)
        # df = df[df["SVTYPE"] != "BND"]
    
        # ✅ SVLEN yoksa INFO'dan çıkar, sonra sayısala çevir
        if "SVLEN" not in df or df["SVLEN"].isna().all():
            df["SVLEN"] = df["INFO"].astype(str).str.extract(r'(?:^|;)SVLEN=(-?\d+)', expand=False)
        df["SVLEN"] = pd.to_numeric(df["SVLEN"], errors="coerce")
    
        # ✅ QUAL sayısal
        df["QUAL"] = pd.to_numeric(df.get("QUAL"), errors="coerce")
    
        # Kaynak etiketle
        df["__SOURCE__"] = tag
    
        # SVTYPE olmayanları at (x ekseninde yanlışları temizler)
        df = df.dropna(subset=["SVTYPE"])
    
        return df
        
    def _grouped_basic_section(title, pairs):
        """
        pairs: list of tuples (vcf_path, 'LABEL'), örn:
          [ (tp_base_path, 'TP-BASE'), (tp_comp_path, 'TP-COMP') ]
          [ (fp_path, 'FP'), (fn_path, 'FN') ]
        """
        dfs = []
        for p, lab in pairs:
            d = _load_and_tag(p, lab)
            if d is not None and not d.empty:
                dfs.append(d)
        if not dfs:
            return html.Div(f"⚠️ {title}: no data.")
    
        data = pd.concat(dfs, ignore_index=True)
    
        # 1) SVTYPE dağılımı: grouped bar
        sv_df = (
            data.dropna(subset=["SVTYPE"])
                .groupby(["__SOURCE__", "SVTYPE"])
                .size()
                .reset_index(name="Count")
        )
        fig_svtype = px.bar(
            sv_df, x="SVTYPE", y="Count", color="__SOURCE__",  color_discrete_sequence=px.colors.qualitative.Vivid,barmode="group",
            title=f"{title} — SVTYPE Distribution (Grouped)"
        )
        fig_svtype.update_layout(margin=dict(l=20, r=20, t=60, b=30), height=420)
    
        # 2) Kromozom dağılımı: grouped bar (chr1..22, X, Y ile sınırlayalım)
        chrom_order = [str(i) for i in range(1,23)] + ["X","Y"]
        chrom_df = data[data["CHROM"].isin(chrom_order)].copy()
        chrom_df["CHROM"] = pd.Categorical(chrom_df["CHROM"], categories=chrom_order, ordered=True)
        chrom_df = (
            chrom_df.groupby(["__SOURCE__", "CHROM"])
                    .size()
                    .reset_index(name="Variant Count")
                    .sort_values(["CHROM", "__SOURCE__"])
        )
        fig_chr = px.bar(
            chrom_df, x="CHROM", y="Variant Count", color="__SOURCE__", color_discrete_sequence=px.colors.qualitative.Vivid,barmode="group",
            title=f"{title} — Chromosome-wise Distribution (Grouped)"
        )
        fig_chr.update_layout(margin=dict(l=20, r=20, t=60, b=30), height=420)
    
        # 3) SVLEN dağılımı: histogram (kaynak renkle), violin (yan yana)
        svlen_df = data.dropna(subset=["SVTYPE", "SVLEN"]).copy()
        # histogram (log gösterim istersen log_y veya x ekseninde transform ekleyebiliriz)
        fig_svlen = px.histogram(
            svlen_df, x="SVLEN", color="__SOURCE__", color_discrete_sequence=px.colors.qualitative.Vivid,nbins=40,
            title=f"{title} — SVLEN Histogram (by Source)"
        )
        fig_svlen.update_layout(margin=dict(l=20, r=20, t=60, b=30), height=420)
    
        # violin: SVTYPE’a göre, source yan yana
        violin_df = svlen_df.copy()
        # log kullanmak istersen:
        # violin_df["logSVLEN"] = np.log10(violin_df["SVLEN"].abs() + 1)
        fig_violin = px.violin(
            violin_df, x="SVTYPE", y="SVLEN", color="__SOURCE__", color_discrete_sequence=px.colors.qualitative.Vivid,box=True, points="all",
            title=f"{title} — SVLEN by SVTYPE (Grouped Violin)"
        )
        fig_violin.update_layout(margin=dict(l=20, r=20, t=60, b=30), height=480)


        all_types = sorted(sv_df["SVTYPE"].unique().tolist())
        preferred = ["DEL", "INS", "DUP", "INV", "BND", "CNV"]
        ordered_types = [t for t in preferred if t in all_types] + [t for t in all_types if t not in preferred]
        
        fig_spider = go.Figure()
        if ordered_types:  # veri varsa çiz
            # Kaynaklara (TP-BASE, TP-COMP, FP, FN) göre pivotla
            pivot = (
                sv_df.pivot_table(index="SVTYPE", columns="__SOURCE__", values="Count", aggfunc="sum")
                    .reindex(index=ordered_types)
                    .fillna(0)
                    .astype(int)
            )
            thetas = ordered_types + [ordered_types[0]]  # poligonu kapatmak için ilkini sona ekle
            palette = px.colors.qualitative.Vivid
        
            for idx, lab in enumerate(sv_df["__SOURCE__"].drop_duplicates().tolist()):
                r_vals = pivot[lab].tolist() if lab in pivot.columns else [0] * len(ordered_types)
                r_closed = r_vals + [r_vals[0]]
                fig_spider.add_trace(go.Scatterpolar(
                    r=r_closed,
                    theta=thetas,
                    fill='toself',
                    name=lab,
                    line=dict(width=2, color=palette[idx % len(palette)]),
                    opacity=0.6
                ))
        
        fig_spider.update_layout(
            polar=dict(radialaxis=dict(visible=True)),
            title=f"{title} — Variant Type Radar Chart",
            showlegend=True,
            margin=dict(l=20, r=20, t=60, b=30),
            height=480
                )
        
        return html.Div([
            html.H4(title, style={'marginTop': '15px'}),
            dcc.Graph(figure=fig_svtype),
            dcc.Graph(figure=fig_spider),   # ← EK
            dcc.Graph(figure=fig_chr),
            dcc.Graph(figure=fig_svlen),
            dcc.Graph(figure=fig_violin)
        ])



    def advanced_truvari_section(title, pairs):
        import plotly.express as px 
        import dash_html_components as html
        import os
        dfs = []
        for p, lab in pairs:
            d = _load_and_tag(p, lab)
            if d is not None and not d.empty:
                dfs.append(d)
        if not dfs:
            return html.Div(f"⚠️ {title}: no data.")
        
        data = pd.concat(dfs, ignore_index=True)
        # --- SANKEY: TP-BASE ➝ SVTYPE ➝ TP-COMP ---
        sankey_df = data.dropna(subset=["CHROM", "SVTYPE", "__SOURCE__"]).copy()
        sankey_df["CHROM"] = "chr" + sankey_df["CHROM"].astype(str)
        
        # CHROM etiketlerini BASE ve COMP'e göre ayır
        sankey_df["CHROM_SRC"] = sankey_df.apply(
            lambda row: f"{row['CHROM']} (BASE)" if row["__SOURCE__"] == "TP-BASE" else f"{row['CHROM']} (COMP)", axis=1
        )
        
        # BASE için: CHROM_SRC → SVTYPE
        base_df = sankey_df[sankey_df["__SOURCE__"] == "TP-BASE"].copy()
        base_df_grouped = base_df.groupby(["CHROM_SRC", "SVTYPE"]).size().reset_index(name="Count")
        
        # COMP için: SVTYPE → CHROM_SRC
        comp_df = sankey_df[sankey_df["__SOURCE__"] == "TP-COMP"].copy()
        comp_df_grouped = comp_df.groupby(["SVTYPE", "CHROM_SRC"]).size().reset_index(name="Count")
        
        # Etiketler
        base_chroms = sorted(base_df_grouped["CHROM_SRC"].unique())
        svtypes = sorted(set(base_df_grouped["SVTYPE"]).union(set(comp_df_grouped["SVTYPE"])))
        comp_chroms = sorted(comp_df_grouped["CHROM_SRC"].unique())
        all_labels = base_chroms + svtypes + comp_chroms
        label_map = {label: i for i, label in enumerate(all_labels)}
        
        # BASE → SVTYPE bağlantıları
        base_source = base_df_grouped["CHROM_SRC"].map(label_map)
        base_target = base_df_grouped["SVTYPE"].map(label_map)
        base_value = base_df_grouped["Count"]
        
        # SVTYPE → COMP bağlantıları
        comp_source = comp_df_grouped["SVTYPE"].map(label_map)
        comp_target = comp_df_grouped["CHROM_SRC"].map(label_map)
        comp_value = comp_df_grouped["Count"]
        
        # SVTYPE'a özel renkler
        svtype_color_map = {sv: px.colors.qualitative.Vivid[i % 10] for i, sv in enumerate(svtypes)}
        base_color = base_df_grouped["SVTYPE"].map(svtype_color_map)
        comp_color = comp_df_grouped["SVTYPE"].map(svtype_color_map)
        
        # Sankey oluştur
        sankey_fig = go.Figure(go.Sankey(
            arrangement="snap",
            node=dict(
                pad=15,
                thickness=20,
                line=dict(color="black", width=0.5),
                label=all_labels,
                color=["#DC143C"] * len(all_labels)
            ),
            link=dict(
                source=pd.concat([base_source, comp_source]),
                target=pd.concat([base_target, comp_target]),
                value=pd.concat([base_value, comp_value]),
                color=pd.concat([base_color, comp_color])
            )
        ))
        
        sankey_fig.update_layout(
            title=f"{title} — TP-BASE → SVTYPE → TP-COMP Sankey",
            font_size=12,
            height=600,
            margin=dict(l=20, r=20, t=60, b=30)
        )

        # ---------- CLUSTERGRAM ----------
        clust_df = data.copy()
        
        # Sadece TP-BASE / TP-COMP kullan (FP/FN vs. gelirse karışmasın)
        clust_df = clust_df[clust_df["__SOURCE__"].isin(["TP-BASE", "TP-COMP"])].copy()
        
        # CHROM "chr" prefiksi (görsel tutarlılık)
        clust_df["CHROM"] = clust_df["CHROM"].astype(str).str.strip()
        clust_df["CHROM"] = clust_df["CHROM"].apply(lambda x: x if x.startswith("chr") else f"chr{x}")
        
        # Çoklu sütun: (__SOURCE__, SVTYPE)
        pivot = clust_df.pivot_table(
            index="CHROM",
            columns=["__SOURCE__", "SVTYPE"],
            aggfunc="size",
            fill_value=0
        )
        
        # Kromozom sıralaması (chr1..22, chrX, chrY) + sadece mevcut olanları al
        chrom_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
        pivot = pivot.reindex([c for c in chrom_order if c in pivot.index])
        
        # Bazı SVTYPE/SOURCE kombinasyonları yoksa kolonları tamlamak isteyebilirsin (opsiyonel):
        # desired_sv = ["DEL", "INS", "DUP", "INV", "BND", "CNV"]
        # desired_sources = ["TP-BASE", "TP-COMP"]
        # import itertools
        # full_cols = pd.MultiIndex.from_product([desired_sources, desired_sv], names=["__SOURCE__", "SVTYPE"])
        # pivot = pivot.reindex(columns=full_cols, fill_value=0)
        
        # Kolon etiketlerini okunur yap: "TP-BASE:DEL" gibi düz string'e çevir
        if isinstance(pivot.columns, pd.MultiIndex):
            nice_cols = [f"{src}:{sv}" for (src, sv) in pivot.columns]
        else:
            # Çok nadir bir durumda MultiIndex düşerse
            nice_cols = pivot.columns.astype(str).tolist()
        
        clustergram_fig = Clustergram(
            data=pivot.values,
            row_labels=pivot.index.tolist(),     # chr1..chrY
            column_labels=nice_cols,             # "TP-BASE:DEL", "TP-COMP:DEL", ...
            color_map="RdBu",
            height=800,
            width=1000,
            display_ratio=[0.85, 0.15],
            # hidden_labels=[],
            color_threshold={"row": 0.5, "col": 0.5}
        )
        
        # --- Tüm x/y eksenlerinde tick label'ları aç ---
        layout_dict = clustergram_fig.layout.to_plotly_json()
        axis_updates = {}
        for key in layout_dict.keys():
            if str(key).startswith("xaxis") or str(key).startswith("yaxis"):
                axis_updates[key] = dict(showticklabels=True)
        
        clustergram_fig.update_layout(
            **axis_updates,
            margin=dict(l=110, r=60, t=50, b=110)  # etiketlere alan
        )
        
    # ------------------ CIRCOS (TP-BASE & TP-COMP) ENTEGRASYONU --------------------
        # TP-BASE / TP-COMP VCF yollarını 'pairs'ten çek
        # --- Genome length dictionaries ---
        # --- Genome lengths ---
        hg19_lengths = {
            "chr1": 249250621, "chr2": 243199373, "chr3": 198022430,
            "chr4": 191154276, "chr5": 180915260, "chr6": 171115067,
            "chr7": 159138663, "chr8": 146364022, "chr9": 141213431,
            "chr10": 135534747, "chr11": 135006516, "chr12": 133851895,
            "chr13": 115169878, "chr14": 107349540, "chr15": 102531392,
            "chr16": 90354753,  "chr17": 81195210,  "chr18": 78077248,
            "chr19": 59128983,  "chr20": 63025520,  "chr21": 48129895,
            "chr22": 51304566,  "chrX": 155270560,  "chrY": 59373566
        }
        hg38_lengths = {
            "chr1": 248956422, "chr2": 242193529, "chr3": 198295559,
            "chr4": 190214555, "chr5": 181538259, "chr6": 170805979,
            "chr7": 159345973, "chr8": 145138636, "chr9": 138394717,
            "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
            "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
            "chr16": 90338345,  "chr17": 83257441,  "chr18": 80373285,
            "chr19": 58617616,  "chr20": 64444167,  "chr21": 46709983,
            "chr22": 50818468,  "chrX": 156040895,  "chrY": 57227415
        }

        # --- Chromosome colors and order ---
        chrom_colors = [
            "#f28e2b", "#4e79a7", "#e15759", "#76b7b2", "#59a14f", "#edc948",
            "#b07aa1", "#ff9da7", "#9c755f", "#bab0ab", "#8cd17d", "#b6992d",
            "#499894", "#d37295", "#fabfd2", "#d4a6c8", "#9f9f9f", "#bcbd22",
            "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#c5b0d5", "#c49c94"
        ]
        chrom_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
        chrom_color_map = dict(zip(chrom_order, chrom_colors))
                       
        tp_base_vcf = next((p for p, lab in pairs if lab == "TP-BASE"), None)
        tp_comp_vcf = next((p for p, lab in pairs if lab == "TP-COMP"), None)
        
        circos_div = html.Div("⚠️ Circos: TP-BASE / TP-COMP not found.")
        if tp_base_vcf and tp_comp_vcf and os.path.exists(tp_base_vcf) and os.path.exists(tp_comp_vcf):
            import tempfile, uuid
            tmp_json = os.path.join(tempfile.gettempdir(), f"circos_truvari_{uuid.uuid4().hex}.json")
            
            # JSON üret
            vcf_to_circos_json_truvari(
                tp_base_vcf, tp_comp_vcf,
                tmp_json,
                svtype_color_map=svtype_color_map,
                chrom_color_map=chrom_color_map,
                hg38_lengths=hg38_lengths,
                hg19_lengths=hg19_lengths,
                bin_size=10_000_000
            )
            
            with open(tmp_json) as f:
                circos_data = json.load(f)
            
            layout = circos_data.get("layout", [])
            layout_ids = {str(b.get("id")) for b in layout}
            all_tracks = circos_data.get("tracks", [])
            
            # Her kombinasyona (örn. "DEL (TP-BASE)") farklı renk
            big_palette = (
                px.colors.qualitative.Vivid
                + px.colors.qualitative.Set2
                + px.colors.qualitative.Plotly
                + px.colors.qualitative.Safe
                + px.colors.qualitative.Bold
                + px.colors.qualitative.Pastel
            )
            def combo_color(name: str) -> str:
                return big_palette[abs(hash(name)) % len(big_palette)]
            
            # BASE ve COMP için ayrı bin listeleri
            base_bins, comp_bins = [], []
            
            for t in all_tracks:
                if str(t.get("type", "")).lower() != "histogram":
                    continue
                name = t.get("name", "")
                lname = name.lower()
                is_base = ("tp-base" in lname) or name.endswith(" (TP-BASE)")
                is_comp = ("tp-comp" in lname) or name.endswith(" (TP-COMP)")
                if not (is_base or is_comp):
                    continue
            
                col = combo_color(name)
                for item in t.get("data", []):
                    blk = str(item.get("block_id"))
                    if blk not in layout_ids:
                        continue
                    try:
                        s = int(item.get("start", 0))
                        e = int(item.get("end", 0))
                        v = int(item.get("value", 0))
                    except Exception:
                        continue
                    if e <= s or v < 0:
                        continue
                    datum = {"block_id": blk, "start": s, "end": e, "value": v, "color": col}
                    (base_bins if is_base else comp_bins).append(datum)
            
            # Yeni, farklı renkler tanımlayalım
            base_color = "#3366ff" # Mavi
            comp_color = "#ff6666" # Kırmızı
            
            tracks = []
            # ✅ TP-BASE = içteki bant (ideogramın içinde)
            if base_bins:
                tracks.append({
                    "type": "HISTOGRAM",
                    "name": "TP-BASE",
                    "data": base_bins,
                    "innerRadius": 160,
                    "outerRadius": 195,
                    "config": {
                        "innerRadius": 160,
                        "outerRadius": 195,
                        "opacity": 1.0,
                        "strokeWidth": 0,
                        "color": base_color # 💡 Renk eklendi
                    }
                })
            
            # ✅ TP-COMP = onun üstünde ikinci iç bant
            if comp_bins:
                tracks.append({
                    "type": "HISTOGRAM",
                    "name": "TP-COMP",
                    "data": comp_bins,
                    "innerRadius": 200,
                    "outerRadius": 235,
                    "config": {
                        "innerRadius": 200,
                        "outerRadius": 235,
                        "opacity": 1.0,
                        "strokeWidth": 0,
                        "color": comp_color # 💡 Renk eklendi
                    }
                })
            
            if not tracks:
                circos_div = html.Div("⚠️ Circos: no histogram data available.")
            else:
                # Kombinasyon legend’i
                combo_names = sorted({
                    t.get("name", "") for t in all_tracks
                    if str(t.get("type", "")).lower() == "histogram"
                })
                legend = html.Div([
                    html.Div([
                        html.Span(style={
                            'backgroundColor': combo_color(lbl),
                            'display': 'inline-block', 'width': '14px', 'height': '14px',
                            'marginRight': '6px', 'borderRadius': '3px'
                        }),
                        html.Span(lbl, style={'marginRight': '12px'})
                    ], style={'display': 'inline-block', 'marginBottom': '4px'})
                    for lbl in combo_names
                ], style={'padding': '6px 0'})
            
                circos_div = html.Div([
                    html.Hr(),
                    html.H4("Circos — TP-BASE & TP-COMP"),
                    legend,
                    dashbio.Circos(
                        id="circos-truvari",
                        layout=layout,
                        tracks=tracks,
                        config={
                            "innerRadius": 260,
                            "outerRadius": 320,
                            "ticks": {"display": True, "spacing": 10_000_000, "labelSpacing": 5},
                            "labelLayout": {"spacing": 60, "radialOffset": 90}
                        },
                        size=800
                    )
                ])
        # -------------------------------------------------------------------------------

        return html.Div([
            html.H4(f"{title} — Advanced Graphs", style={'marginTop': '15px'}),
            dcc.Graph(figure=sankey_fig),
            dcc.Graph(figure=clustergram_fig),
            circos_div,  # Circos block
        ])
    # ✅ TP-BASE vs TP-COMP (yan yana)
    tp_block_basic = _grouped_basic_section(
        "TP (TP-BASE vs TP-COMP)",
        [(tru_paths.get("tp_base"), "TP-BASE"), (tru_paths.get("tp_comp"), "TP-COMP")]
    )
    # ✅ Advanced TP grafikleri (sadece TP için)
    tp_block_advanced = advanced_truvari_section(
        "TP (TP-BASE vs TP-COMP)",
        [(tru_paths.get("tp_base"), "TP-BASE"), (tru_paths.get("tp_comp"), "TP-COMP")]
    )    
    # ❌ FP vs FN (yan yana)
    err_block = _grouped_basic_section(
        "FP vs FN",
        [(tru_paths.get("fp"), "FP"), (tru_paths.get("fn"), "FN")]
    )

    visuals_block = html.Div([
        tp_block_basic,
        html.Hr(),
        tp_block_advanced,
        html.Hr(),
        err_block  # FP/FN basic only
    ])    
    
    # Dosya listesi istenmiyor -> boş
    return status_div, html.Ul([]), metrics_preview, visuals_block
