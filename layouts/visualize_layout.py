from dash import html, dcc, dash_table
import dash_bootstrap_components as dbc

def get_visualize_layout():
    return dbc.Container([
        html.Div([
            html.Div(id="visual_output"),
            dcc.Store(id='circos-data-store'),
            dcc.Store(id='svtype-colors-store'),
            dcc.Store(id='manhattan-plot-type-store'),   # ← tracks which Manhattan mode is active
            html.Div(id="circos-container"),
            html.Div(id="default-circos-output"),
            dcc.Store(id='manhattan-store'),

            html.Div(
                id="manhattan-controls-container",
                style={"display": "none"},
                children=[

                    # ── QUAL-based Manhattan controls ──────────────────────
                    html.Div(
                        id="manhattan-qual-controls",
                        children=[
                            html.Div([
                                html.Label(
                                    "Significance Threshold (−log₁₀(P)):",
                                    style={'marginTop': '10px', 'fontWeight': 'bold'}
                                ),
                                # Info tooltip
                                html.Span(
                                    " ⓘ",
                                    id="qual-threshold-info",
                                    style={
                                        'cursor': 'pointer',
                                        'color': '#0984e3',
                                        'fontSize': '14px'
                                    }
                                ),
                                dbc.Tooltip(
                                    [
                                        html.P(
                                            "The threshold is applied to −log₁₀(P), where P is derived from the "
                                            "VCF QUAL score using the formula P = 10^(−QUAL / 200).",
                                            style={'marginBottom': '4px'}
                                        ),
                                        html.P(
                                            "Default value: 6. Variants above this threshold are highlighted "
                                            "in red as 'Points of Interest'.",
                                            style={'marginBottom': '4px'}
                                        ),
                                        html.P(
                                            "Increasing the threshold makes the filter more stringent; "
                                            "only variants with higher QUAL scores will be highlighted."
                                        ),
                                    ],
                                    target="qual-threshold-info",
                                    placement="right",
                                    style={'maxWidth': '350px'}
                                ),
                            ], style={'display': 'flex', 'alignItems': 'center', 'gap': '4px'}),

                            dcc.Slider(
                                id='manhattan-threshold-slider',
                                min=0, max=10, step=0.5, value=6,
                                marks={i: str(i) for i in range(0, 11)},
                                tooltip={"placement": "bottom", "always_visible": True}
                            ),

                            # Dynamic label showing what current threshold means
                            html.Div(
                                id='manhattan-qual-threshold-label',
                                style={
                                    'fontSize': '12px',
                                    'color': '#555',
                                    'marginTop': '4px',
                                    'fontStyle': 'italic'
                                }
                            ),

                            html.Div([
                                html.Span("📌 ", style={'fontSize': '13px'}),
                                html.Span(
                                    "Variants above the dashed red line are highlighted as 'Points of Interest'.",
                                    style={'fontSize': '12px', 'color': '#555'}
                                )
                            ], style={'marginTop': '6px', 'marginBottom': '10px'}),
                        ]
                    ),

                    # ── SVLEN-based Manhattan controls ─────────────────────
                    html.Div(
                        id="manhattan-svlen-controls",
                        style={"display": "none"},
                        children=[
                            html.Div([
                                html.Label(
                                    "SV Length Threshold (log₁₀ scale):",
                                    style={'marginTop': '10px', 'fontWeight': 'bold'}
                                ),
                                html.Span(
                                    " ⓘ",
                                    id="svlen-threshold-info",
                                    style={
                                        'cursor': 'pointer',
                                        'color': '#0984e3',
                                        'fontSize': '14px'
                                    }
                                ),
                                dbc.Tooltip(
                                    [
                                        html.P(
                                            "The threshold is applied in log₁₀(SVLEN) space, consistent with "
                                            "the logarithmic Y-axis of the plot.",
                                            style={'marginBottom': '4px'}
                                        ),
                                        html.P(
                                            "Default value: 4 → highlights SVs with SVLEN > 10,000 bp.",
                                            style={'marginBottom': '4px'}
                                        ),
                                        html.Ul([
                                            html.Li("Slider value 3 → SVLEN > 1,000 bp"),
                                            html.Li("Slider value 4 → SVLEN > 10,000 bp  (default)"),
                                            html.Li("Slider value 5 → SVLEN > 100,000 bp"),
                                            html.Li("Slider value 6 → SVLEN > 1,000,000 bp"),
                                        ], style={'fontSize': '12px', 'paddingLeft': '16px'}),
                                    ],
                                    target="svlen-threshold-info",
                                    placement="right",
                                    style={'maxWidth': '380px'}
                                ),
                            ], style={'display': 'flex', 'alignItems': 'center', 'gap': '4px'}),

                            dcc.Slider(
                                id='manhattan-svlen-threshold-slider',
                                min=1, max=7, step=0.5, value=4,
                                marks={
                                    1: '10 bp',
                                    2: '100 bp',
                                    3: '1 kb',
                                    4: '10 kb',
                                    5: '100 kb',
                                    6: '1 Mb',
                                    7: '10 Mb',
                                },
                                tooltip={"placement": "bottom", "always_visible": True}
                            ),

                            # Dynamic label
                            html.Div(
                                id='manhattan-svlen-threshold-label',
                                style={
                                    'fontSize': '12px',
                                    'color': '#555',
                                    'marginTop': '4px',
                                    'fontStyle': 'italic'
                                }
                            ),

                            html.Div([
                                html.Span("📌 ", style={'fontSize': '13px'}),
                                html.Span(
                                    "Variants above the dashed red line are highlighted as 'Points of Interest'.",
                                    style={'fontSize': '12px', 'color': '#555'}
                                )
                            ], style={'marginTop': '6px', 'marginBottom': '10px'}),
                        ]
                    ),

                    # ── SVTYPE filter (shared by both plot modes) ──────────
                    html.Br(),
                    html.Label(
                        "Filter by SV Type:",
                        style={'marginTop': '10px', 'fontWeight': 'bold'}
                    ),
                    html.Div(
                        "Select one or more SV types to display. Deselecting all types will clear the plot.",
                        style={'fontSize': '12px', 'color': '#888', 'marginBottom': '4px'}
                    ),
                    dcc.Dropdown(
                        id='manhattan-svtype-selector',
                        options=[],
                        value=[],
                        multi=True,
                        placeholder="Select one or more SV types",
                        style={"width": "300px"}
                    ),
                ]
            ),

            html.Div(id="default-dashbio-manhattanplot")
        ])
    ], fluid=True)
def get_manhattan_controls(slider_id, svtype_id, prefix=""):
    return html.Div([
        html.Div([
            html.Label(
                "Significance Threshold (−log₁₀(P)):",
                style={'marginTop': '10px', 'fontWeight': 'bold'}
            ),
            html.Span(
                " ⓘ",
                id=f"{prefix}qual-threshold-info",
                style={'cursor': 'pointer', 'color': '#0984e3', 'fontSize': '14px'}
            ),
            dbc.Tooltip(
                [
                    html.P(
                        "The threshold is applied to −log₁₀(P), where P is derived from the "
                        "VCF QUAL score using the formula P = 10^(−QUAL / 200).",
                        style={'marginBottom': '4px'}
                    ),
                    html.P(
                        "Default value: 6. Variants above this threshold are highlighted "
                        "in red as 'Points of Interest'.",
                        style={'marginBottom': '4px'}
                    ),
                    html.P(
                        "Increasing the threshold makes the filter more stringent; "
                        "only variants with higher QUAL scores will be highlighted."
                    ),
                ],
                target=f"{prefix}qual-threshold-info",
                placement="right",
                style={'maxWidth': '350px'}
            ),
        ], style={'display': 'flex', 'alignItems': 'center', 'gap': '4px'}),

        dcc.Slider(
            id=slider_id,
            min=0, max=10, step=0.5, value=6,
            marks={i: str(i) for i in range(0, 11)},
            tooltip={"placement": "bottom", "always_visible": True}
        ),

        html.Div(
            id=f"{prefix}qual-threshold-label",
            style={'fontSize': '12px', 'color': '#555', 'marginTop': '4px', 'fontStyle': 'italic'}
        ),

        html.Div([
            html.Span("📌 ", style={'fontSize': '13px'}),
            html.Span(
                "Variants above the dashed red line are highlighted as 'Points of Interest'.",
                style={'fontSize': '12px', 'color': '#555'}
            )
        ], style={'marginTop': '6px', 'marginBottom': '10px'}),

        html.Label(
            "Filter by SV Type:",
            style={'marginTop': '10px', 'fontWeight': 'bold'}
        ),
        html.Div(
            "Select one or more SV types to display. Deselecting all types will clear the plot.",
            style={'fontSize': '12px', 'color': '#888', 'marginBottom': '4px'}
        ),
        dcc.Dropdown(
            id=svtype_id,
            options=[],
            value=[],
            multi=True,
            placeholder="Select one or more SV types",
            style={"width": "300px"}
        ),
    ])