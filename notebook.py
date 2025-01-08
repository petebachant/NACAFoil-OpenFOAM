import marimo

__generated_with = "0.9.33"
app = marimo.App(width="medium")


@app.cell
def _():
    import pandas as pd
    return (pd,)


@app.cell
def _(pd):
    df = pd.read_csv("processed/all-simulated.csv").sort_values("alpha_deg")
    df
    return (df,)


@app.cell
def _(df):
    df.plot(x="alpha_deg", y="cl", backend="plotly")
    return


@app.cell
def _(df):
    df["cl_cd"] = df["cl"] / df["cd"]
    df.plot(x="alpha_deg", y="cl_cd", backend="plotly")
    return


@app.cell
def __():
    return


if __name__ == "__main__":
    app.run()
