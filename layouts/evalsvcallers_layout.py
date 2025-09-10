from dash import html, dcc
import dash_bootstrap_components as dbc

def get_evalsvcallers_layout():
    shared_input_style = {
        'fontFamily': '"Times New Roman", Times, serif',
        'width': '200px',
        'marginLeft': '0px'
    }

    shared_label_style = {
        'fontFamily': '"Times New Roman", Times, serif',
        'width': '250px',
        'display': 'inline-block'
    }

    shared_row_style = {
        'marginBottom': '8px'
    }

    return html.Div([
  #  html.H3("EvalSVcallers", style={'padding': '1rem', 'fontSize': '24px', 'fontFamily': '"Times New Roman", Times, serif', 'fontWeight': 'bold'}),

     html.H3("Upload Caller VCF File", style={'marginTop': '20px', 'fontFamily': '"Times New Roman", Times, serif', 'fontSize': '20px','whiteSpace': 'nowrap'}),
    # Hidden conversion section
    html.Div(id='conversion-section', children=[
    
        # Upload raw VCFs
        dcc.Upload(
            id='eval-upload',
            children=html.Div(['📎 Drag and Drop or Select File']),
            style={
                'width': '100%', 'height': '60px', 'lineHeight': '60px',
                'borderWidth': '1px', 'borderStyle': 'dashed', 'borderRadius': '5px',
                'textAlign': 'center', 'margin': '0px', 'fontFamily': '"Times New Roman", Times, serif'
            },
            multiple=False
        ),
    
        #html.Div(id='uploaded-files-list', style={'marginTop': '20px'}),
    #    html.Br(),
    
        html.H3("Select Caller Tool for Conversion", style={'marginTop': '20px','fontFamily': '"Times New Roman", Times, serif', 'fontSize': '20px','whiteSpace': 'nowrap'}),
    
        dcc.Dropdown(
            id='caller_tool',
            options=[{"label": tool, "value": tool} for tool in [
                "BASIL-ANISE", "BICseq2", "BreakDancer", "BreakSeek", "BreakSeq2", "Breakway",
                "CNVnator", "Control-FREEC", "CREST", "DELLY", "DIGTYPER", "DINUMT","ERDS",
                "FermiKit", "forestSV", "GASVpro","GenomeStrip", "GRIDSS","HGT-ID", "Hydra-sv", "iCopyDAV",
                "inGAP-sv", "ITIS","laSV","Lumpy", "Manta", "MATCHCLIP", "Meerkat","MELT", "metaSV", "MindTheGap",
                "Mobster","OncoSNP-Seq", "Pindel", "PBHoney-NGM","PBHoney","PennCNV-Seq", "PRISM", "RAPTR", "readDepth","RetroSeq",
                "Sniffles", "Socrates", "SoftSearch", "SoftSV", "SoloDel", "SvABA",
                "SVDetect", "SVfinder", "SVseq2", "SVIM","Tangram","Tea","TEMP","TIDDIT", "Ulysses", "Vaquita","VariationHunter",
                "VirusFinder2", "Wham"
            ]],
            placeholder="Select a caller tool",
            style={'fontFamily': '"Times New Roman", Times, serif', 'width': '100%'}
        ),
    
        dbc.Button(
            'Convert VCF File',
            id='convert-button',
            n_clicks=0,
            color="primary",
            className="mt-2",
            style={'fontFamily': '"Times New Roman", Times, serif'}
        )
    
       # html.Div(id='convert-status', style={'marginTop': '10px', 'fontFamily': '"Times New Roman", Times, serif'}),
    #    html.Div(id='converted-file-list', style={'marginTop': '10px', 'fontFamily': '"Times New Roman", Times, serif'}),
    ]), 


    html.Div(id='converted-upload-section', children=[
   #     html.H3("Upload Converted VCF File", style={
   #         'marginTop': '0px', 'fontFamily': '"Times New Roman", Times, serif',
   #         'fontSize': '20px', 'whiteSpace': 'nowrap'
    #    }),
        dcc.Upload(
            id='upload-converted-data',
            children=html.Div(['📎 Drag and Drop or Select File']),
            style={
                'width': '400px', 'height': '60px', 'lineHeight': '60px',
                'borderWidth': '1px', 'borderStyle': 'dashed', 'borderRadius': '5px',
                'textAlign': 'center', 'margin': '10px',
                'fontFamily': '"Times New Roman", Times, serif'
            },
            multiple=False
        )
    ], style={'display': 'none'}),

  #  html.Div(id='converted-files-list', style={'marginTop': '10px', 'fontFamily': '"Times New Roman", Times, serif'}),
    html.Br(),

    # Reference selection
    html.H3("Select Reference Type", style={'fontFamily': '"Times New Roman", Times, serif', 'fontSize': '20px','whiteSpace': 'nowrap'}),
    dcc.RadioItems(
        id="reference_choice",
        options=[
            {"label": "Embedded NA12878", "value": "N"},
            {"label": "Embedded Sim-A", "value": "A"},
            {"label": "Univar SV Katalog", "value": "univar"},
            {"label": "Upload Custom Reference", "value": "custom"}
        ],
        value="N",
        inline=True,
        labelStyle={'display': 'block', 'fontFamily': '"Times New Roman", Times, serif'}
    ),


    html.Div(id='reference_upload_section', children=[

        dcc.Upload(
            id='dynamic-upload',
            children=html.Div(['📎 Drag and Drop or Select File']),
            style={
                'width': '100%', 'height': '60px', 'lineHeight': '60px',
                'borderWidth': '1px', 'borderStyle': 'dashed', 'borderRadius': '5px',
                'textAlign': 'center', 'margin': '10px', 'fontFamily': '"Times New Roman", Times, serif'
            },
            multiple=False
        )
    ], style={'width': '100%'}),


        
 #   html.Div(id="reference_upload_section", style={"marginTop": "10px"}),
    html.Br(),

    # Parameter mode toggle
    html.H3("Parameter Selection Mode", style={'fontFamily': '"Times New Roman", Times, serif', 'fontSize': '20px','whiteSpace': 'nowrap'}),
    dcc.RadioItems(
        id='param_mode',
        options=[{'label': 'Basic', 'value': 'basic'}, {'label': 'Advanced', 'value': 'advanced'}],
        value='basic',
        inline=True,
        labelStyle={'display': 'inline-block', 'marginRight': '10px', 'fontFamily': '"Times New Roman", Times, serif'}
    ),

    # Basic Parameters
    html.Div(id='basic_params', children=[
        html.Div([
            html.Label('SV Type (-st):', style=shared_label_style),
            dcc.Input(id='param_st_basic', type='text', placeholder='default: ALL e.g. DEL', style=shared_input_style)
        ], style=shared_row_style),

        html.Div([
            html.Label('Minimum SV Length (-l):', style=shared_label_style),
            dcc.Input(id='param_l_basic', type='number', min=0, placeholder='default: 50', style=shared_input_style)
        ], style=shared_row_style),

        html.Div([
            html.Label('Maximum SV Length (-xl):', style=shared_label_style),
            dcc.Input(id='param_xl_basic', type='number', min=0, placeholder='default: 2000000', style=shared_input_style)
        ], style=shared_row_style),

        html.Div([
            html.Label('Min REF Len (-rl):', style=shared_label_style),
            dcc.Input(id='param_rl_basic', type='number', min=0, placeholder='default: 30', style=shared_input_style)
        ], style=shared_row_style),

        html.Div([
            html.Label('Max REF Len (-rxl):', style=shared_label_style),
            dcc.Input(id='param_rxl_basic', type='number', min=0, placeholder='default: 2000000', style=shared_input_style)
        ], style=shared_row_style),

        html.Div([
            html.Label('Minimum Overlap (-mo):', style=shared_label_style),
            dcc.Input(id='param_mo_basic', type='number', min=0, max=1, step=0.01, placeholder='default: 0.5', style=shared_input_style)
        ], style=shared_row_style),

        html.Div([
            html.Label('Output TF (-of): (0: no output, 1: TP calls, 2: FP calls, 3: TP and FP calls)', style=shared_label_style),
            dcc.Input(id='param_of_basic', type='number', min=0, max=3, placeholder='default: 3', style=shared_input_style)
        ], style=shared_row_style),
    ], style={'marginBottom': '20px','display': 'block'}),

    # Advanced Parameters
    html.Div(id='advanced_params', children=[
 
        html.Div([
            html.Label('Target Chromosome (-c):', style=shared_label_style),
            dcc.Input(id='param_c', type='text', placeholder='default: all', style=shared_input_style)
        ], style=shared_row_style),

        html.Div([
            html.Label('Minimum Reads (-mr):', style=shared_label_style),
            dcc.Input(id='param_mr', type='number', min=0, placeholder='default: 0', style=shared_input_style)
        ], style=shared_row_style),

        html.Div([
            html.Label('Read Set (-rs):', style=shared_label_style),
            dcc.Input(id='param_rs', type='number', min=1, max=2, placeholder='default: 1', style=shared_input_style)
        ], style=shared_row_style),

        html.Div([
            html.Label('Minimum Insertion (-mins):', style=shared_label_style),
            dcc.Input(id='param_mins', type='number', min=0, placeholder='default: 200', style=shared_input_style)
        ], style=shared_row_style),
        
        html.Div([
            html.Label('Evaluate Genotype (-eg): [default: false]', style=shared_label_style),
            dcc.Checklist(id='param_eg', options=[{'label': ' True', 'value': 'True'}], style=shared_label_style)
        ], style=shared_row_style),

        html.Div([
            html.Label('Evaluate Breakpoints (-eb): [default: false]', style=shared_label_style),
            dcc.Checklist(id='param_eb', options=[{'label': ' True', 'value': 'True'}], style=shared_label_style)
        ], style=shared_row_style),

        html.Div([
            html.Label('Insertion Eval (-i): [default: false]', style=shared_label_style),
            dcc.Checklist(id='param_i', options=[{'label': ' True', 'value': 'True'}], style=shared_label_style)
        ], style=shared_row_style),

        html.Div([
            html.Label('Include Y Chromosome (-y): [default: false]', style=shared_label_style),
            dcc.Checklist(id='param_y', options=[{'label': ' True', 'value': 'True'}], style=shared_label_style)
        ], style=shared_row_style),

        html.Div([
            html.Label('Subtract ME (-sm): [default: false]', style=shared_label_style),
            dcc.Checklist(id='param_sm', options=[{'label': ' True', 'value': 'True'}], style=shared_label_style)
        ], style=shared_row_style),

        html.Div([
            html.Div([
                html.Label('Parent 1 (-p1):', style={**shared_label_style, 'flex': '0 0 150px'}),
                dcc.Upload(id='param_parent1', children=html.Button('Upload Parent 1'), multiple=False, style={'fontFamily': '"Times New Roman", Times, serif','width': '200px'})
            ], style={'display': 'flex', 'alignItems': 'center'}),
        ], style=shared_row_style),

        html.Div([
            html.Div([
                html.Label('Parent 2 (-p2):', style={**shared_label_style, 'flex': '0 0 150px'}),
                dcc.Upload(id='param_parent2', children=html.Button('Upload Parent 2'), multiple=False, style={'fontFamily': '"Times New Roman", Times, serif','width': '200px'})
            ], style={'display': 'flex', 'alignItems': 'center'}),
        ], style=shared_row_style),

        html.Div([
            html.Div([
                html.Label('Region BED (-rb):', style={**shared_label_style, 'flex': '0 0 150px'}),
                dcc.Upload(id='param_rb', children=html.Button('Upload BED'), multiple=False, style={'fontFamily': '"Times New Roman", Times, serif','width': '200px'})
            ], style={'display': 'flex', 'alignItems': 'center'}),
        ], style=shared_row_style)
    ]), 
   # html.Br(),
    dbc.Button('Run EvalSVcallers', id='eval-button', n_clicks=0, color="primary", className="mt-2", style={'fontFamily': '"Times New Roman", Times, serif'}),
    #html.Div(id='eval-status', style={'marginTop': '15px', 'fontFamily': '"Times New Roman", Times, serif'}),
   # dcc.Loading(
    #    id="loading-eval-wrapper",
       # type="default",
      #  children=html.Div([
        #    html.Div(id='eval-status', style={'marginTop': '20px'}),
         #   html.Div(id='eval-output-list', style={'marginTop': '10px', 'fontFamily': '"Times New Roman", Times, serif'})
        #])
   # )
    #html.Div(id='eval-output-list', style={'marginTop': '10px', 'fontFamily': '"Times New Roman", Times, serif'}),
    #html.Div(id='eval-output', style={'marginTop': '20px', 'fontFamily': '"Times New Roman", Times, serif'})
])