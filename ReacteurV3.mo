package ReacteurV3
  import SI = Modelica.SIunits;

  model reacteurV3 "Ideal flow source that produces a prescribed mass flow with prescribed specific enthalpy, mass fraction and trace substances"
    import Modelica.Media.Interfaces.Choices.IndependentVariables;
    extends PartialFlowSource;
    parameter Medium.Temperature T = Medium.T_default "Fixed value of temperature";
    parameter Medium.MassFraction X[Medium.nX] = Medium.X_default "Fixed value of composition";
    //parameter Medium.ExtraProperty C[Medium.nC](quantity = Medium.extraPropertiesNames) = fill(0, Medium.nC) "Fixed values of trace substances" annotation(Evaluate = true, Dialog(enable = not use_C_in and Medium.nC > 0));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort;
    Medium.MassFlowRate m_flow "Fixed mass flow rate going out of the fluid port";
    Medium.SaturationProperties sat "State vector to compute saturation properties";
    Medium.SpecificEnthalpy h_l = Medium.bubbleEnthalpy(sat);
    Medium.SpecificEnthalpy h_v = Medium.dewEnthalpy(sat) "specific enthalpy of vapour";
    //SI.EnthalpyFlowRate Hb_flow "Enthalpy flow across boundaries or energy source/sink";
    SI.HeatFlowRate Qb_flow "Heat flow across boundaries or energy source/sink";
  equation
    Modelica.Fluid.Utilities.checkBoundary(Medium.mediumName, Medium.substanceNames, Medium.singleState, true, X, "MassFlowSource_h");
    medium.T = T;
    sat.Tsat = medium.T;
    sat.psat = Medium.saturationPressure(medium.T);
    heatPort.T = medium.T;
    m_flow = -Qb_flow / (h_v - h_l);
    port.m_flow = m_flow;
    Qb_flow = heatPort.Q_flow;
    medium.Xi = X[1:Medium.nXi];
    //port.C_outflow = C;
    port.h_outflow = h_v;
  end reacteurV3;

  partial model PartialFlowSource "Partial component source with one fluid connector"
    import Modelica.Constants;
    //    parameter Integer nPorts = 0 "Number of ports" annotation(Dialog(connectorSizing = true));
    replaceable package Medium = Modelica.Media.Interfaces.PartialMedium "Medium model within the source" annotation(choicesAllMatching = true);
    Medium.BaseProperties medium "Medium in the source";
    Modelica.Fluid.Interfaces.FluidPort_b port(redeclare package Medium = Medium, m_flow(max = 0, min = -Constants.inf)) annotation(Placement(transformation(extent = {{90, 10}, {110, -10}})));
  equation
    port.p = medium.p;
    port.Xi_outflow = medium.Xi;
  end PartialFlowSource;
end ReacteurV3;