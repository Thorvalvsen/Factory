model TurbineSolaireV2
  import SI = Modelica.SIunits;

  model Turbine008
    inner Modelica.Fluid.System system(p_ambient = 1.2e5, T_ambient = 300, allowFlowReversal = false);
    replaceable package MediumL = Modelica.Media.Water.StandardWater;
    replaceable package MediumV = Modelica.Media.Water.WaterIF97OnePhase_ph;
    final constant Integer nbEtage = 200;
    ReacteurV3.reacteurV3 reacteur(T = 380, redeclare package Medium = MediumL);
    baseClasses.Concentrateur concentrateur(surface = 3);
    baseClasses.TableSolaire tablesolaire;
    Modelica.Fluid.Sources.FixedBoundary boundary(redeclare package Medium = MediumV, T = 380, p = 100000, nPorts = 1);
    turboMachine.turbineAxiale Etage[nbEtage](diamMobile = 0.008, diamHub = 0.007, angleGuide = 50, angleMobile = 5, nu = 0.9, redeclare package Medium = MediumV);
    //turboMachine.turbineAxiale Etage2(diamMobile = 0.008, diamHub = 0.007, angleGuide = 50, angleMobile = 5, nu = 0.9, redeclare package Medium = MediumV);
    Entrainement.multiplexeur multi(nbEntree = nbEtage);
    Entrainement.recuperateur Generateur(N = 10000);
  equation
    connect(tablesolaire.solar.y, concentrateur.valeursoleil);
    connect(concentrateur.heatPort, reacteur.heatPort);
    connect(reacteur.port, Etage[1].port_a);
    for i in 1:nbEtage - 1 loop
      connect(Etage[i].port_b, Etage[i + 1].port_a);
      connect(Etage[i].shaft, multi.arbreE[i]);
    end for;
    connect(Etage[nbEtage].port_b, boundary.ports[1]);
    connect(Etage[nbEtage].shaft, multi.arbreE[nbEtage]);
    connect(multi.arbreS, Generateur.shaft);
  end Turbine008;

  model Turbine2
    inner Modelica.Fluid.System system(p_ambient = 1.2e5, T_ambient = 300, allowFlowReversal = false);
    replaceable package MediumL = Modelica.Media.Water.StandardWater;
    replaceable package MediumV = Modelica.Media.Water.IdealSteam;
    ReacteurV3.reacteurV3 reacteur(T = 400, redeclare package Medium = MediumL);
    baseClasses.Concentrateur concentrateur(surface = 2);
    baseClasses.TableSolaire tablesolaire;
    Modelica.Fluid.Sources.FixedBoundary boundary(redeclare package Medium = MediumV, T = 380, p = 100000, nPorts = 1);
    turboMachine.turbineAxiale Etage1(diamMobile = 0.008, diamHub = 0.007, angleGuide = 50, angleMobile = 5, nu = 0.9, redeclare package Medium = MediumV);
    Entrainement.recuperateur Moteur1(N = 13000);
  equation
    connect(tablesolaire.solar.y, concentrateur.valeursoleil);
    connect(concentrateur.heatPort, reacteur.heatPort);
    connect(reacteur.port, Etage1.port_a);
    connect(Etage1.port_b, boundary.ports[1]);
    connect(Etage1.shaft, Moteur1.shaft);
  end Turbine2;

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

  package Entrainement
    model EntrainementCV "Entrainement avec couple fonction de la vitesse"
      extends entrainement;
      parameter SI.AngularVelocity omegaS "Vitesse de synchronisme";
      parameter SI.Torque Tmax "Couple Maxi";
      Real g "Glissement";
      Real gmax = 0.925 "Glissement a couple maxi";
    equation
      g = (omega - omegaS) / omegaS;
      Tn = 2 * Tmax / gmax * g / (1 + g ^ 2 / gmax ^ 2);
      shaft.tau = Tn;
    end EntrainementCV;

    model EntrainementVC "Entrainement a Vitesse constante"
      import SI = Modelica.SIunits;
      extends entrainement;
      parameter SI.AngularVelocity N = 1500 "Rotation en tr/min";
    equation
      omega = N * 2 * Modelica.Constants.pi / 60;
    end EntrainementVC;

    model EntrainementVPC "Entrainement a Vitesse et puissance constante"
      Modelica.Mechanics.Rotational.Interfaces.Flange_a shaft;
      parameter SI.AngularVelocity N "Rotation en tr/min";
      parameter SI.Power Pmeca;
      SI.Angle phi(start = 0.001) "Angle de rotation de l'arbre de la turbine";
    equation
      phi = shaft.phi;
      der(phi) = SI.Conversions.from_rpm(N);
      shaft.tau = -Pmeca / SI.Conversions.from_rpm(N);
    end EntrainementVPC;

    model recuperateur "Recuperateur d'energie"
      extends entrainement;
      SI.Energy Wm "Energie mecanique";
      parameter SI.AngularVelocity N = 1500 "Rotation en tr/min";
    equation
      omega = N * 2 * Modelica.Constants.pi / 60;
      Wm = phi * shaft.tau;
    end recuperateur;

    model multiplexeur
      parameter Integer nbEntree = 1;
      SI.Angle phi "Angle de rotation";
      SI.Torque tau;
      Modelica.Mechanics.Rotational.Interfaces.Flange_a arbreE[nbEntree];
      Modelica.Mechanics.Rotational.Interfaces.Flange_b arbreS;
    equation
      phi = arbreS.phi;
      sum(arbreE.tau) = tau;
      for i in 1:nbEntree loop
        arbreE[i].phi = phi;
      end for;
      arbreS.tau = tau;
    end multiplexeur;

    model spacer
      SI.Angle phi "Angle de rotation";
      SI.Torque tau;
      Modelica.Mechanics.Rotational.Interfaces.Flange_a arbreE;
      Modelica.Mechanics.Rotational.Interfaces.Flange_b arbreS;
    equation
      phi = arbreS.phi;
      arbreE.tau = tau;
      arbreE.phi = phi;
      arbreS.tau = tau;
    end spacer;

    partial model entrainement "Moteur avec transmission de puissance"
      import SI = Modelica.SIunits;
      SI.Angle phi(start = 0.001) "Angle de rotation de l'arbre de la turbine";
      SI.AngularVelocity omega "vitesse de rotation de l'arbre";
      SI.Power Pm "Puissance mecanique transmise";
      Modelica.Mechanics.Rotational.Interfaces.Flange_a shaft;
    equation
      phi = shaft.phi;
      omega = der(phi);
      Pm = omega * shaft.tau;
    end entrainement;
  end Entrainement;

  package baseClasses
    connector HeatPort "Connecteur pour transmettre la chaleur recuperee"
      extends Modelica.Thermal.HeatTransfer.Interfaces.HeatPort;
    end HeatPort;

    model TableSolaire
      Modelica.Blocks.Sources.TimeTable solar(table = [0, 0; 3600, 100; 7200, 150; 10800, 200; 14400, 350; 18000, 350; 21600, 200; 25200, 150; 28800, 100]);
    end TableSolaire;

    model Concentrateur "Parabole avec une surface parametrable et ensoleillement variable"
      import baseClasses.HeatPort;
      baseClasses.HeatPort heatPort;
      parameter Boolean use_soleil = true;
      parameter Modelica.SIunits.Area surface = 1 "Surface de la parabole";
      Modelica.Blocks.Interfaces.RealInput valeursoleil(unit = "W/m2");
    equation
      heatPort.Q_flow = -surface * valeursoleil;
    end Concentrateur;

    model generateur
      parameter SI.AngularVelocity N = 1500 "Rotation en tr/min";
      // parameter Integer nb_port;
      SI.Angle phi "Angle de rotation de l'arbre de la turbine";
      SI.AngularVelocity omega "vitesse de rotation de l'arbre";
      SI.Power Pm "Puissance mecanique transmise";
      Modelica.Mechanics.Rotational.Interfaces.Flange_a shaft;
      SI.Power sumPm;
    equation
      omega = N * 2 * Modelica.Constants.pi / 60;
      //for i in 1:nb_port loop
      phi = shaft.phi;
      omega = der(phi);
      Pm = omega * shaft.tau;
      sumPm = sumPm + Pm;
      //end for;
    end generateur;
  end baseClasses;
end TurbineSolaireV2;