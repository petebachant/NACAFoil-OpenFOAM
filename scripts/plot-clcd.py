"""Plot lift-to-drag ratio."""

import pandas as pd
import plotly.graph_objects as go

fig = go.Figure()

df_exp = pd.read_csv("processed/NACA0012_6e6_Ladson_180grit.csv")
df_exp["cl_cd"] = df_exp["cl"] / df_exp["cd"]
df_exp = df_exp[df_exp.alpha_deg > 0]
df_sim = pd.read_csv("processed/all-simulated.csv").sort_values("alpha_deg")
df_sim["cl_cd"] = df_sim.cl / df_sim.cd
fig.add_trace(
    go.Scatter(
        x=df_sim.alpha_deg, y=df_sim["cl_cd"], mode="lines+markers", name="CFD"
    )
)
fig.add_trace(
    go.Scatter(
        x=df_exp.alpha_deg, y=df_exp["cl_cd"], mode="markers", name="Exp"
    )
)
fig.update_layout(
    xaxis=dict(title="Angle of attack (deg)"),
    yaxis=dict(title="Lift-to-drag ratio"),
    showlegend=True,
    title=None,
    margin=dict(t=50)
)
fig.write_json("figures/naca0012-clcd.json")
