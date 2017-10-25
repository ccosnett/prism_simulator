(* Mathematica Package *)
(* Created by Mathematica Plugin for IntelliJ IDEA *)

(* :Title: prism *)
(* :Context: prism` *)
(* :Author: johncosnett *)
(* :Date: 2017-10-25 *)

(* :Package Version: 0.1 *)
(* :Mathematica Version: *)
(* :Copyright: (c) 2017 johncosnett *)
(* :Keywords: *)
(* :Discussion: simple ray tracing and internal reflections using snell's law and the law of reflection*)


BeginPackage["prism`"(*,{"snell1`","snell2`","botPoint1`","lefPoint1`","righPoint1`","botPoint2`","lefPoint2`","righPoint2`"}*)];
(** main construct **)
s::usage = "s[]";
glass::usage = "glass[\[Theta],n]";
internalRay0::usage = "internalRay0[\[Theta],n]";

sourceRay::usage = "sourceRay[\[Theta], n]";
sourceRay2::usage = "sourceRay2[\[Theta],n]";
sourceRay3::usage = "sourceRay3[\[Theta],n]";

ray1::usage = "ray1[\[Theta],n]";
ray2::usage = "ray2[\[Theta],n]";
ray3::usage = "ray3[\[Theta],n]";
ray4::usage = "ray4[\[Theta],n]";
ray5::usage = "ray5[\[Theta],n]";
ray6::usage = "ray6[\[Theta],n]";

txtThetaT2::usage = "txtThetaT2[\[Theta],n]";
txtThetaI1::usage = "txtThetaI1[\[Theta],n]";
txtN::usage = "txtN[\[Theta],n]";
txtDelta::usage = "txtDelta[\[Theta],n]";



(** tools and metric functions **)
delta::usage = "delta[\[Theta],n]";
r2D::usage = "r2D[]";


(** summoning points **)

rightPointOne::usage = "rightPointOne[\[Theta],n]";
rightPointTwo::usage = "rightPointTwo[\[Theta],n]";
bottomPointOne::usage = "bottomPointOne[\[Theta],n]";
bottomPointTwo::usage = "bottomPointTwo[\[Theta],n]";
leftPointOne::usage = "leftPointOne[\[Theta],n]";
leftPointTwo::usage = "leftPointTwo[\[Theta],n]";

Begin["`Private`"];


txtDelta[\[Theta]_, n_] :=
 g[Text["\[Delta] = " <> t[AccountingForm[delta[\[Theta], n], 3]] <>
    "\[Degree]", {1.01, 1}, Left]]
txtN[\[Theta]_, n_] := g[Text["n = " <> t[AccountingForm[N[n,3], 3]], {1.01, 1.3}, Left]]
txtThetaI1[\[Theta]_, n_] :=
 g[Text["\!\(\*SubscriptBox[\(\[Theta]\), \(i1\)]\) = " <>
    t[AccountingForm[r2D[\[Theta]], 3]] <> "\[Degree]", {1.01, 1.6},
   Left]]
txtThetaT2[\[Theta]_, n_] :=
 g[Text["\!\(\*SubscriptBox[\(\[Theta]\), \(t2\)]\) = " <>
    t[AccountingForm[thetaOuterRight[\[Theta], n], 3]] <>
    "\[Degree]", {1.01, 1.9}, Left]]

delta[\[Theta]_, n_] := thetaOuterLeft[\[Theta], n] + thetaOuterRight[\[Theta], n] - r2D[\[Pi]/3]


s=Show[base, rightInterface, leftInterface,  (*strikerLeft,*)##,rainj]&;



ray1[\[Theta]_, n_] :=
  g[{gg, l[{rightInterfacePoint[\[Theta], n][[1, 1]],
      bottomPointOne[\[Theta], n][[1, 1]]}]}];
ray2[\[Theta]_, n_] :=
  g[{gg, l[{bottomPointOne[\[Theta], n][[1, 1]],
      leftPointOne[\[Theta], n][[1, 1]]}]}];
ray3[\[Theta]_, n_] :=
  g[{gg, l[{leftPointOne[\[Theta], n][[1, 1]],
      rightPointOne[\[Theta], n][[1, 1]]}]}];
ray4[\[Theta]_, n_] :=
  g[{gg, l[{rightPointOne[\[Theta], n][[1, 1]],
      bottomPointTwo[\[Theta], n][[1, 1]]}]}];
ray5[\[Theta]_, n_] :=
  g[{gg, l[{bottomPointTwo[\[Theta], n][[1, 1]],
      leftPointTwo[\[Theta], n][[1, 1]]}]}];
ray6[\[Theta]_, n_] :=
  g[{gg, l[{leftPointTwo[\[Theta], n][[1, 1]],
      rightPointTwo[\[Theta], n][[1, 1]]}]}];

thetaOuterLeft[\[Theta]_, n_] :=
  Module[{m1, m2, x1, y1, x2, y2, point1, point2},
   point1 = sourcePoint[\[Theta]][[1, 1]];
   point2 = strikerLeftCoordinate;
   x1 = point1[[1]];
   y1 = point1[[2]];

   x2 = point2[[1]];
   y2 = point2[[2]];
   m1 = -(1/Sqrt[3]);
   m2 = (y2 - y1)/(x2 - x1);
   r2D[ArcTan[-((m1 - m2)/(1 + m1 m2))] // N]
   ];
thetaOuterRight[\[Theta]_, n_] :=
  Module[{m1, m2, x1, y1, x2, y2, point1, point2},
   point1 = sourcePoint2[\[Theta], n][[1, 1]];
   point2 = rightInterfacePoint[\[Theta], n][[1, 1]];
   x1 = point1[[1]];
   y1 = point1[[2]];

   x2 = point2[[1]];
   y2 = point2[[2]];
   m1 = 1/Sqrt[3];
   m2 = (y2 - y1)/(x2 - x1);
   r2D[ArcTan[(m1 - m2)/(1 + m1 m2)] // N]
   ];


(*rayLast*)
(************************************************)
(************************************************)
rightPointTwo[\[Theta]_, n_] := Module[
  {xy, x1, y1, m, x, y},
  xy = leftPointTwo[\[Theta], n][[1, 1]];
  x1 = xy[[1]];
  y1 = xy[[2]];
  m = mLefMirror2[\[Theta], n];

  x = (m x1 + Sqrt[3] - y1)/(Sqrt[3] + m);
  y = (m - m x1 + y1)/(1 + m/Sqrt[3]);
  g[p[{x, y}]]
  ];




mLefMirror2[\[Theta]_, n_] := Tan[((11 Pi)/6) + thetaBottomMirror2NormalDifference[\[Theta], n]];
thetaBottomMirror2NormalDifference[\[Theta]_, n_] := Module[{m1, m2},
  m1 = -(1/Sqrt[3]);
  m2 = mBotMirror2[\[Theta], n];
  ArcTan[(m1 - m2)/(1 + m1 m2)]
  ];

















leftPointTwo[\[Theta]_, n_] := Module[
  {xy, x1, y1, m, x, y},
  xy = bottomPointTwo[\[Theta], n][[1,1]];
  x1 = xy[[1]];
  y1 = xy[[2]];
  m = mBotMirror2[\[Theta], n];


  x = (-m x1 - Sqrt[3] + y1)/(Sqrt[3] - m);
  y = (y1 - m - m x1)/(1 - m/Sqrt[3]);

  g[p[{x, y}]]
  ];

mBotMirror2[\[Theta]_,n_]:=-mRighMirror2[\[Theta],n];


bottomPointTwo[\[Theta]_, n_] := Module[

  {xy, x1, y1, m, x, y},

  xy = rightPointOne[\[Theta], n][[1, 1]];
  x1 = xy[[1]];
  y1 = xy[[2]];
  m = mRighMirror2[\[Theta], n];

  x = x1 - y1/m;
  y = 0;
  g[p[{x, y}]]
  ];

mRighMirror2[\[Theta]_, n_] := Tan[\[Pi]/6 + thetaLefMirror1NormalDifference[\[Theta], n]];

thetaLefMirror1NormalDifference[\[Theta]_, n_] := Module[{m1, m2},
  m1 = 1/Sqrt[3];
  m2 = mLefMirror1[\[Theta], n];
  ArcTan[(m1 - m2)/(1 + m1 m2)]
  ];




rightPointOne[\[Theta]_, n_] := Module[
  {xy, x1, y1, m, x, y},
  xy = leftPointOne[\[Theta], n][[1, 1]];
  x1 = xy[[1]];
  y1 = xy[[2]];
  m = mLefMirror1[\[Theta], n];

  x = (m x1 + Sqrt[3] - y1)/(Sqrt[3] + m);
  y = (m - m x1 + y1)/(1 + m/Sqrt[3]);
  g[p[{x, y}]]
  ];




mLefMirror1[\[Theta]_, n_] := Tan[((11 Pi)/6) + thetaBottomMirror1NormalDifference[\[Theta], n]];
thetaBottomMirror1NormalDifference[\[Theta]_, n_] := Module[{m1, m2},
  m1 = -(1/Sqrt[3]);
  m2 = mBotMirror1[\[Theta], n];
  ArcTan[(m1 - m2)/(1 + m1 m2)]
  ];
leftPointOne[\[Theta]_, n_] := Module[
  {xy, x1, y1, m, x, y},
  xy = bottomPointOne[\[Theta], n][[1,1]];
  x1 = xy[[1]];
  y1 = xy[[2]];
  m = mBotMirror1[\[Theta], n];


  x = (-m x1 - Sqrt[3] + y1)/(Sqrt[3] - m);
  y = (y1 - m - m x1)/(1 - m/Sqrt[3]);

  g[p[{x, y}]]
  ];

mBotMirror1[\[Theta]_,n_]:=-mRighMirror1[\[Theta],n];

bottomPointOne[\[Theta]_, n_] := Module[

  {xy, x1, y1, m, x, y},

  xy = rightInterfacePoint[\[Theta], n][[1, 1]];
  x1 = xy[[1]];
  y1 = xy[[2]];
  m = mRighMirror1[\[Theta], n];

  x = x1 - y1/m;
  y = 0;
  g[p[{x, y}]]
  ];

mRighMirror1[\[Theta]_, n_] := Tan[\[Pi]/6 + thetaInternalNormalDifference[\[Theta], n]];
thetaInternalNormalDifference[\[Theta]_, n_] := Module[{m1, m2},
  m1 = 1/Sqrt[3];
  m2 = mInternal[\[Theta], n];
  ArcTan[(m1 - m2)/(1 + m1 m2)]
  ];





mInternal[\[Theta]_, n_] := Tan[\[Theta]t1[\[Theta], n] - \[Pi]/6];
\[Theta]t1[\[Theta]_, n_] := ArcSin[Sin[\[Theta]]/n];

p=Point;
g=Graphics;
l=Line;

rightInterfacePoint[\[Theta]_, n_] := g[p[{xRightSource[\[Theta], n], yRightSource[\[Theta], n]}]];
xRightSource[\[Theta]_, n_] := (Sqrt[3] - y1Source + mInternal[\[Theta], n] x1Source)/(mInternal[\[Theta], n] + Sqrt[3]) ;
yRightSource[\[Theta]_, n_] := (mInternal[\[Theta], n] - mInternal[\[Theta], n] x1Source + y1Source)/(1 + mInternal[\[Theta], n]/Sqrt[3]);
y1Source=Sqrt[3]/3;
x1Source=-2/3;
mInternal[\[Theta]_, n_] := Tan[\[Theta]t1[\[Theta], n] - \[Pi]/6];
\[Theta]t1[\[Theta]_, n_] := ArcSin[Sin[\[Theta]]/n];

(************************************************)
(************************************************)


(*snell2*)
(************************************************)
(************************************************)
sourceRay3[\[Theta]_, n_] := g[{Darker[Green], l[{rightPointTwo[\[Theta], n][[1, 1]], sourcePoint3[\[Theta], n][[1, 1]]}]}];
sourcePoint3[\[Theta]_, n_] := g[Point[{Cos[-thetaT3[\[Theta], n] + \[Pi]/6], Sin[-thetaT3[\[Theta], n] + \[Pi]/6]} + rightPointTwo[\[Theta], n][[1, 1]]]];
thetaT3[\[Theta]_, n_] := Module[{a}, a = ArcSin[n*Sin[rightInterfaceAngle2[\[Theta], n]]]; If[Internal`RealValuedNumericQ[a], a, Pi/2]];
rightInterfaceAngle2[\[Theta]_, n_] := (*r2D@*)((Pi/2)-ArcTan[(mLefMirror2[\[Theta], n] - (-Sqrt[3]))/(1 + mLefMirror2[\[Theta], n] (-Sqrt[3]))]);



(*mInternal -> mLefMirror2*)
(*rightInterfacePoint -> rightPointTwo*)
(*rightInterfaceAngle -> rightInterfaceAngle2*)
(*thetaT2 -> thetaT3*)
(*sourcePoint2 -> sourcePoint3*)
(*sourceRay2 -> sourceRay3*)



rightPointTwo[\[Theta]_, n_] := Module[
  {xy, x1, y1, m, x, y},
  xy = leftPointTwo[\[Theta], n][[1, 1]];
  x1 = xy[[1]];
  y1 = xy[[2]];
  m = mLefMirror2[\[Theta], n];

  x = (m x1 + Sqrt[3] - y1)/(Sqrt[3] + m);
  y = (m - m x1 + y1)/(1 + m/Sqrt[3]);
  g[p[{x, y}]]
  ];
mLefMirror2[\[Theta]_, n_] := Tan[((11 Pi)/6) + thetaBottomMirror2NormalDifference[\[Theta], n]];
thetaBottomMirror2NormalDifference[\[Theta]_, n_] := Module[{m1, m2},
  m1 = -(1/Sqrt[3]);
  m2 = mBotMirror2[\[Theta], n];
  ArcTan[(m1 - m2)/(1 + m1 m2)]
  ];

leftPointTwo[\[Theta]_, n_] := Module[
  {xy, x1, y1, m, x, y},
  xy = bottomPointTwo[\[Theta], n][[1,1]];
  x1 = xy[[1]];
  y1 = xy[[2]];
  m = mBotMirror2[\[Theta], n];


  x = (-m x1 - Sqrt[3] + y1)/(Sqrt[3] - m);
  y = (y1 - m - m x1)/(1 - m/Sqrt[3]);

  g[p[{x, y}]]
  ];

mBotMirror2[\[Theta]_,n_]:=-mRighMirror2[\[Theta],n];


bottomPointTwo[\[Theta]_, n_] := Module[

  {xy, x1, y1, m, x, y},

  xy = rightPointOne[\[Theta], n][[1, 1]];
  x1 = xy[[1]];
  y1 = xy[[2]];
  m = mRighMirror2[\[Theta], n];

  x = x1 - y1/m;
  y = 0;
  g[p[{x, y}]]
  ];

mRighMirror2[\[Theta]_, n_] := Tan[\[Pi]/6 + thetaLefMirror1NormalDifference[\[Theta], n]];

thetaLefMirror1NormalDifference[\[Theta]_, n_] := Module[{m1, m2},
  m1 = 1/Sqrt[3];
  m2 = mLefMirror1[\[Theta], n];
  ArcTan[(m1 - m2)/(1 + m1 m2)]
  ];




rightPointOne[\[Theta]_, n_] := Module[
  {xy, x1, y1, m, x, y},
  xy = leftPointOne[\[Theta], n][[1, 1]];
  x1 = xy[[1]];
  y1 = xy[[2]];
  m = mLefMirror1[\[Theta], n];

  x = (m x1 + Sqrt[3] - y1)/(Sqrt[3] + m);
  y = (m - m x1 + y1)/(1 + m/Sqrt[3]);
  g[p[{x, y}]]
  ];




mLefMirror1[\[Theta]_, n_] := Tan[((11 Pi)/6) + thetaBottomMirror1NormalDifference[\[Theta], n]];
thetaBottomMirror1NormalDifference[\[Theta]_, n_] := Module[{m1, m2},
  m1 = -(1/Sqrt[3]);
  m2 = mBotMirror1[\[Theta], n];
  ArcTan[(m1 - m2)/(1 + m1 m2)]
  ];
leftPointOne[\[Theta]_, n_] := Module[
  {xy, x1, y1, m, x, y},
  xy = bottomPointOne[\[Theta], n][[1,1]];
  x1 = xy[[1]];
  y1 = xy[[2]];
  m = mBotMirror1[\[Theta], n];


  x = (-m x1 - Sqrt[3] + y1)/(Sqrt[3] - m);
  y = (y1 - m - m x1)/(1 - m/Sqrt[3]);

  g[p[{x, y}]]
  ];

mBotMirror1[\[Theta]_,n_]:=-mRighMirror1[\[Theta],n];

bottomPointOne[\[Theta]_, n_] := Module[

  {xy, x1, y1, m, x, y},

  xy = rightInterfacePoint[\[Theta], n][[1, 1]];
  x1 = xy[[1]];
  y1 = xy[[2]];
  m = mRighMirror1[\[Theta], n];

  x = x1 - y1/m;
  y = 0;
  g[p[{x, y}]]
  ];

mRighMirror1[\[Theta]_, n_] := Tan[\[Pi]/6 + thetaInternalNormalDifference[\[Theta], n]];
thetaInternalNormalDifference[\[Theta]_, n_] := Module[{m1, m2},
  m1 = 1/Sqrt[3];
  m2 = mInternal[\[Theta], n];
  ArcTan[(m1 - m2)/(1 + m1 m2)]
  ];





mInternal[\[Theta]_, n_] := Tan[\[Theta]t1[\[Theta], n] - \[Pi]/6];
\[Theta]t1[\[Theta]_, n_] := ArcSin[Sin[\[Theta]]/n];

p=Point;
g=Graphics;
l=Line;

rightInterfacePoint[\[Theta]_, n_] := g[p[{xRightSource[\[Theta], n], yRightSource[\[Theta], n]}]];
xRightSource[\[Theta]_, n_] := (Sqrt[3] - y1Source + mInternal[\[Theta], n] x1Source)/(mInternal[\[Theta], n] + Sqrt[3]) ;
yRightSource[\[Theta]_, n_] := (mInternal[\[Theta], n] - mInternal[\[Theta], n] x1Source + y1Source)/(1 + mInternal[\[Theta], n]/Sqrt[3]);
y1Source=Sqrt[3]/3;
x1Source=-2/3;
mInternal[\[Theta]_, n_] := Tan[\[Theta]t1[\[Theta], n] - \[Pi]/6];
\[Theta]t1[\[Theta]_, n_] := ArcSin[Sin[\[Theta]]/n];


(************************************************)
(************************************************)


(*snell1*)
(************************************************)
(************************************************)
sourceRay2[\[Theta]_, n_] := g[{Red, l[{rightInterfacePoint[\[Theta], n][[1, 1]], sourcePoint2[\[Theta], n][[1, 1]]}]}]
sourcePoint2[\[Theta]_, n_] := g[Point[{Cos[-thetaT2[\[Theta], n] + \[Pi]/6], Sin[-thetaT2[\[Theta], n] + \[Pi]/6]} + rightInterfacePoint[\[Theta], n][[1, 1]]]]
thetaT2[\[Theta]_, n_] := Module[{a}, a = ArcSin[n*Sin[rightInterfaceAngle[\[Theta], n]]]; If[Internal`RealValuedNumericQ[a], a, Pi/2]];

rightInterfaceAngle[\[Theta]_, n_] := (*r2D@*)((Pi/2)-ArcTan[(mInternal[\[Theta], n] - (-Sqrt[3]))/(1 + mInternal[\[Theta], n] (-Sqrt[3]))]);

rightInterfacePoint[\[Theta]_, n_] := g[p[{xRightSource[\[Theta], n], yRightSource[\[Theta], n]}]];



xRightSource[\[Theta]_, n_] := (Sqrt[3] - y1Source + mInternal[\[Theta], n] x1Source)/(mInternal[\[Theta], n] + Sqrt[3]) ;
yRightSource[\[Theta]_, n_] := (mInternal[\[Theta], n] - mInternal[\[Theta], n] x1Source + y1Source)/(1 + mInternal[\[Theta], n]/Sqrt[3]);
y1Source=Sqrt[3]/3;
x1Source=-2/3;

internalRay0[\[Theta]_, n_] := g[{Red, Line[{ strikerLeftCoordinate, rightInterfacePoint[\[Theta], n][[1, 1]]}]}];
mInternal[\[Theta]_, n_] := Tan[\[Theta]t1[\[Theta], n] - \[Pi]/6]
\[Theta]t1[\[Theta]_, n_] := ArcSin[Sin[\[Theta]]/n]
g = Graphics;
l=Line;
p = Point;
strikerLeftCoordinate = {-2/3, Sqrt[3]/3};
r2D = (180.0/\[Pi]) # &;
d2R = (\[Pi]/180.0) # &;
(************************************************)
(************************************************)


sourceRay[\[Theta]_,n_] := g[{Red, Line[{sourcePoint[\[Theta]][[1, 1]], strikerLeftCoordinate}]}];
mInternal[\[Theta]_, n_] := Tan[\[Theta]t1[\[Theta], n] - \[Pi]/6]
\[Theta]t1[\[Theta]_, n_] := ArcSin[Sin[\[Theta]]/n]
sourcePoint[\[Theta]_] := g[Point[{Cos[\[Theta] + (5 \[Pi])/6], Sin[\[Theta] + (5 \[Pi])/6]} + strikerLeftCoordinate]]
rainj = PlotRange -> {{-2, 2}, {-0.5, 2}};
strikerLeftCoordinate = {-2/3, Sqrt[3]/3};
mythicalStrikerRightCoordinate={2/3,Sqrt[3]/3};
strikerLeft = g[p[{-2/3, Sqrt[3]/3}]];
strikerRight = g[p[{-2/3, Sqrt[3]/3}]];
base = g[{l[{leftCorner, rightCorner}]}];
rightInterface = g[{l[{leftCorner, topCorner}]}];
leftInterface = g[{l[{topCorner, rightCorner}]}];
glass[\[Theta]_,n_] := g[{Darker[White,((n - 1)/2.0)], Polygon[{leftCorner, topCorner, rightCorner}]}];

l=Line;
cornerPoints = g[{PointSize[.02], leftCorner//p, rightCorner//p, topCorner//p}];
g = Graphics;
p = Point;
leftCorner = {-1, 0} ;
rightCorner = {1, 0} ;
topCorner = {0, Sqrt[3]} ;
n=1.7;
t=ToString;
ggg=Darker[Green];
gg=Green;



End[];
EndPackage[];