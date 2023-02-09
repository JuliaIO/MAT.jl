model FallingBodyBox
  inner Modelica.Mechanics.MultiBody.World world annotation(
    Placement(visible = true, transformation(origin = {-70, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Mechanics.MultiBody.Parts.BodyBox bodyBox(r = {1, 0, 0})  annotation(
    Placement(visible = true, transformation(origin = {10, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Mechanics.MultiBody.Joints.FreeMotion freeMotion annotation(
    Placement(visible = true, transformation(origin = {-30, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
  connect(world.frame_b, freeMotion.frame_a) annotation(
    Line(points = {{-60, 10}, {-40, 10}}));
  connect(freeMotion.frame_b, bodyBox.frame_a) annotation(
    Line(points = {{-20, 10}, {0, 10}}, color = {95, 95, 95}));

annotation(
    uses(Modelica(version = "4.0.0")));
end FallingBodyBox;
