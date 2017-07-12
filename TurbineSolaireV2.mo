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
    BancEssai.recuperateur Moteur1(N = 13000);
  equation
    connect(tablesolaire.solar.y, concentrateur.valeursoleil);
    connect(concentrateur.heatPort, reacteur.heatPort);
    connect(reacteur.port, Etage1.port_a);
    connect(Etage1.port_b, boundary.ports[1]);
    connect(Etage1.shaft, Moteur1.shaft);
  end Turbine2;

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