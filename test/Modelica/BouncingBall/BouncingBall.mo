model BouncingBall
  parameter Real eff=0.77 "coefficient of restitution";
  parameter Real grav=9.81 "gravity acceleration";
  Real height(fixed=true, start=111) "height of ball";
  Real vel(fixed=true) "velocity of ball";
  Boolean flying(fixed=true, start=true) "true, if ball is flying";
  Boolean impact;
  Real v_new(fixed=true);
  Integer foo;

equation
  impact = height <= 0.0;
  foo = if impact then 1 else 2;
  der(vel) = if flying then -grav else 0;
  der(height) = vel;

  when {height <= 0.0 and vel <= 0.0, impact} then
  //when {impact} then
    v_new = if edge(impact) then -eff*pre(vel) else 0;
    flying = v_new > 0;
    reinit(vel, v_new);
  end when;

  //copy files from C:/Users/BenConrad/AppData/Local/Temp/OpenModelica/OMEdit/BouncingBall
  annotation(
    experiment(StartTime=0, StopTime=1, Tolerance=1e-3, Interval=0.1)
  );

end BouncingBall;
