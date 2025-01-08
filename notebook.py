import marimo

__generated_with = "0.10.9"
app = marimo.App(width="medium")


@app.cell
def _():
    import pandas as pd
    return (pd,)


@app.cell
def _(pd):
    df = pd.read_csv("processed/all-simulated.csv")
    df
    return (df,)


@app.cell
def _(df):
    df.sort_values("alpha_deg").plot(x="alpha_deg", y="cl", backend="plotly")
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
