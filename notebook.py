import marimo

__generated_with = "0.23.16"
app = marimo.App(width="medium", layout_file="layouts/notebook.grid.json")


@app.cell
def _():
    import marimo as mo
    import pandas as pd
    import plotly.graph_objects as go

    # Resolves to this directory when run with a kernel, and to the served
    # base URL when running in the browser via WebAssembly. Paths beneath it
    # mirror the project layout, so a cell reads the same either way.
    public = mo.notebook_location() / "public"

    slider = mo.ui.slider(0, 20, 1, label="Angle of attack (deg)")
    slider
    return go, mo, pd, public, slider


@app.cell
def _(go, pd, public, slider):
    df = pd.read_csv(
        str(public / "processed" / "all-simulated.csv")
    ).sort_values("alpha_deg")
    df["cl_cd"] = df.cl / df.cd
    dfe = pd.read_csv(
        str(public / "processed" / "NACA0012_6e6_Ladson_180grit.csv")
    )
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
    return


@app.cell
def _(mo, public, slider):
    # Load flow snapshot for this angle of attack
    mo.image(
        str(
            public
            / "figures"
            / f"naca0012-re2e5-aoa-{slider.value}-umag.png"
        )
    )
    return


if __name__ == "__main__":
    app.run()
