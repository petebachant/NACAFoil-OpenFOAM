import marimo

__generated_with = "0.10.9"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    slider = mo.ui.slider(0, 20, 1, label="Angle of attack (deg)")
    slider
    return mo, slider


@app.cell
def _(slider):
    import pandas as pd
    import plotly.graph_objects as go

    df = pd.read_csv("processed/all-simulated.csv").sort_values("alpha_deg")
    df["cl_cd"] = df.cl / df.cd
    dfe = pd.read_csv("processed/NACA0012_6e6_Ladson_180grit.csv")
    dfe["cl_cd"] = dfe.cl / dfe.cd
    dfe = dfe[dfe.alpha_deg > 0]

    fig = go.Figure()

    fig.add_trace(go.Scatter(x=df.alpha_deg, y=df.cl_cd, name="CFD", mode="lines+markers"))
    fig.add_trace(go.Scatter(x=dfe.alpha_deg, y=dfe.cl_cd, name="exp"))
    fig.update_layout(xaxis_title="Angle of attack (deg)", yaxis_title="$C_l/C_d$")
    fig.add_vline(x=slider.value, line_dash="dash")
    fig
    return df, dfe, fig, go, pd


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
