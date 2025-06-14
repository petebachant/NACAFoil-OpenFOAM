import marimo

__generated_with = "0.10.9"
app = marimo.App(width="medium", layout_file="layouts/notebook.grid.json")


@app.cell
def _():
    import marimo as mo
    import plotly.io as pio

    pio.renderers.default = "notebook"

    slider = mo.ui.slider(0, 20, 1, label="Angle of attack (deg)")
    slider
    return mo, slider


@app.cell
def _(slider):
    import pandas as pd
    import plotly.graph_objects as go
    from calkit.datasets import read_dataset

    df = read_dataset("processed/all-simulated.csv").sort_values("alpha_deg")
    df["cl_cd"] = df.cl / df.cd
    dfe = read_dataset("processed/NACA0012_6e6_Ladson_180grit.csv")
    dfe["cl_cd"] = dfe.cl / dfe.cd
    dfe = dfe[dfe.alpha_deg > 0]

    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x=df.alpha_deg, y=df.cl_cd, name="CFD", mode="lines+markers"
        )
    )
    fig.add_trace(go.Scatter(x=dfe.alpha_deg, y=dfe.cl_cd, name="Exp."))
    fig.update_layout(
        xaxis_title="Angle of attack (deg)", yaxis_title="$C_l/C_d$"
    )
    fig.add_vline(x=slider.value, line_dash="dash")
    fig
    return df, dfe, fig, go, pd, read_dataset


@app.cell
def _(mo, slider):
    # Load flow snapshot for this angle of attack
    import calkit

    mo.image(
        calkit.read_file(f"figures/naca0012-re2e5-aoa-{slider.value}-umag.png")
    )
    return (calkit,)


if __name__ == "__main__":
    app.run()
