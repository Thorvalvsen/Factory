package Entrainement
  import SI = Modelica.SIunits;

  model EntrainementCV "Entrainement avec couple fonction de la vitesse"
    extends baseClasse.Entrainement;
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
    extends baseClasse.Entrainement;
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
    extends baseClasse.Entrainement;
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

  package baseClasse
    partial model Entrainement "Moteur avec transmission de puissance"
      import SI = Modelica.SIunits;
      SI.Angle phi(start = 0.001) "Angle de rotation de l'arbre de la turbine";
      SI.AngularVelocity omega "vitesse de rotation de l'arbre";
      SI.Power Pm "Puissance mecanique transmise";
      Modelica.Mechanics.Rotational.Interfaces.Flange_a shaft;
    equation
      phi = shaft.phi;
      omega = der(phi);
      Pm = omega * shaft.tau;
    end Entrainement;
  end baseClasse;
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})));
end Entrainement;