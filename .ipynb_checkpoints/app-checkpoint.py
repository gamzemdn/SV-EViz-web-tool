**app.py**

```python
from dash import Dash, html, dcc, Input, Output
import dash_bootstrap_components as dbc
from dash_bootstrap_templates import load_figure_template

# Import refactored layouts
from survivor_layout import get_survivor_layout
from evalsvcallers_layout import get_evalsvcallers_layout
from visualize import get_visualize_layout
from metrics import get_metrics_layout

load_figure_template('LITERA')
app = Dash(__name__, suppress_callback_exceptions=True, external_stylesheets=[dbc.themes.BOOTSTRAP])

navbar = dbc.Navbar(
    dbc.Container([
        dbc.NavbarBrand("Variant Comparison and Visualization Tool", href="/")
    ]),
    color="primary",
    dark=True
)

app.layout = html.Div([
    navbar,
    dcc.Tabs(id='main-tabs', value='tab-welcome', children=[
        dcc.Tab(label='Welcome', value='tab-welcome'),
        dcc.Tab(label='Comparison', value='tab-comparison'),
        dcc.Tab(label='Visualization', value='tab-visualization'),
        dcc.Tab(label='Metrics', value='tab-metrics'),
    ]),
    html.Div(id='tabs-content')
])

@app.callback(
    Output('tabs-content', 'children'),
    Input('main-tabs', 'value')
)
def render_main_tab(tab):
    if tab == 'tab-welcome':
        return html.Div([
            html.H2("Welcome to the Variant Comparison and Visualization Tool!"),
            html.P("Select an option from the tabs.")
        ])
    elif tab == 'tab-comparison':
        return html.Div([
            html.H4("Select Comparison Type:"),
            dcc.RadioItems(
                id='comparison-type',
                options=[
                    {'label': 'SURVIVOR', 'value': 'survivor'},
                    {'label': 'EvalSVcallers', 'value': 'evalsvcallers'}
                ],
                value='survivor',
                labelStyle={'display': 'block'}
            ),
            html.Div(id='comparison-content')
        ])
    elif tab == 'tab-visualization':
        return get_visualize_layout()
    elif tab == 'tab-metrics':
        return get_metrics_layout()

@app.callback(
    Output('comparison-content', 'children'),
    Input('comparison-type', 'value')
)
def render_comparison_content(selection):
    if selection == 'survivor':
        return get_survivor_layout()
    elif selection == 'evalsvcallers':
        return get_evalsvcallers_layout()

if __name__ == "__main__":
    app.run_server(debug=True, port=8040)
```

