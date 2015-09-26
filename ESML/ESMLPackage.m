(* ::Package:: *)

(************************************************************************)
(* This file was generated automatically by the Mathematica front end.  *)
(* It contains Initialization cells from a Notebook file, which         *)
(* typically will have the same name as this file except ending in      *)
(* ".nb" instead of ".m".                                               *)
(*                                                                      *)
(* This file is intended to be loaded into the Mathematica kernel using *)
(* the package loading commands Get or Needs.  Doing so is equivalent   *)
(* to using the Evaluate Initialization Cells menu command in the front *)
(* end.                                                                 *)
(*                                                                      *)
(* DO NOT EDIT THIS FILE.  This entire file is regenerated              *)
(* automatically each time the parent Notebook file is saved in the     *)
(* Mathematica front end.  Any changes you make to this file will be    *)
(* overwritten.                                                         *)
(************************************************************************)



(* :Title: ESMLPackage *) 
(* :Authors: Vitor V. Vasconcelos & Fernando P. Santos *)
(* :Summary: Summary goes here. *)
(* :Context: ESML`ESMLPackage` *)
(* :Package version: 1.0 *)
(* :History:  Version 1.0 September 16 2015 *)
(* :Mathematica version: 10.2.0 for Mac OS X x86 (64-bit) (July 29, 2015) & 10.2.0 for Microsft Windows (64-bit) (July 7, 2015)*)
(* :Discussion: Give more details here.*)


BeginPackage["ESML`ESMLPackage`"];


(* :Code Section (Call Unprotect and ClearAll): *)


Unprotect["`*"]; 
ClearAll["`*"];


(* :Usage Messages: *)
TransitionMatrixSML::usage="TransitionMatrixSML[s,T,Z] represents the transition matrix between a set of Configurations of Interest (CoI) in the long run of a one-step Markov process at the edges of the phase space (zero order approximation) as well as the set itself. The process is i=(\!\(\*SubscriptBox[\(i\), \(1\)]\),...,\!\(\*SubscriptBox[\(i\), \(s\)]\)) such that \!\(\*UnderoverscriptBox[\(\[Sum]\), \(k = 1\), \(s\)]\)\!\(\*SubscriptBox[\(i\), \(k\)]\)=Z and T[fs,ts,k] is the transition probability from configuration (0,...,\!\(\*SubscriptBox[\(i\), \(fs\)]\)=k,...,\!\(\*SubscriptBox[\(i\), \(ts\)]\)=Z-k,...,0) to (0,...,\!\(\*SubscriptBox[\(i\), \(fs\)]\)=k-1,...,\!\(\*SubscriptBox[\(i\), \(ts\)]\)=Z-k+1,...,0) in each time step. The output is of the form {CoI, TransitionMatrix}, where TransitionMatrix is a sparse array and CoI contains a list {{CoI_1,DDrift_1,diffusion_1}, ...,} where CoI_i is an array {\!\(\*SubscriptBox[\(i\), \(1\)]\),...,\!\(\*SubscriptBox[\(i\), \(s\)]\)} with the position of the CoI, DDrift_i and diffusion_i are nothing but tags.\n"
TransitionMatrixESML::usage="TransitionMatrixESML[s,T,Z] represents the transition matrix between a set of Configurations of Interest (CoI) in the long run of a one-step Markov process at the edges of the phase space (first order approximation) as well as the set itself. The process is i=(\!\(\*SubscriptBox[\(i\), \(1\)]\),...,\!\(\*SubscriptBox[\(i\), \(s\)]\)) such that \!\(\*UnderoverscriptBox[\(\[Sum]\), \(k = 1\), \(s\)]\)\!\(\*SubscriptBox[\(i\), \(k\)]\)=Z and T[fs,ts,k] is the transition probability from configuration (0,...,\!\(\*SubscriptBox[\(i\), \(fs\)]\)=k,...,\!\(\*SubscriptBox[\(i\), \(ts\)]\)=Z-k,...,0) to (0,...,\!\(\*SubscriptBox[\(i\), \(fs\)]\)=k-1,...,\!\(\*SubscriptBox[\(i\), \(ts\)]\)=Z-k+1,...,0) in each time step. The output is of the form {CoI, TransitionMatrix}, where TransitionMatrix is a sparse array and CoI contains a list {{CoI_1,DDrift_1,diffusion_1}, ...,} where CoI_i is an array {\!\(\*SubscriptBox[\(i\), \(1\)]\),...,\!\(\*SubscriptBox[\(i\), \(s\)]\)} with the position of the CoI, DDrift_i and diffusion_i are, respectively, the derivative of the drift (1st Kramers-Moyal term) and the diffusion (2nd Kramers-Moyal term) at the CoI.\n"
StationaryDistributionESML::usage="StationaryDistributionESML[s,Z,ConfigurationsOfInterest,TransitionMatrix] represents an estimation of the stationary distribution of a one-step Markov process at the edges of the phase space. ConfigurationsOfInterest contains a list {{CoI_1,DDrift_1,diffusion_1}, ...,} where CoI_i is an array {\!\(\*SubscriptBox[\(i\), \(1\)]\),...,\!\(\*SubscriptBox[\(i\), \(s\)]\)} with the position of the Configurations of Interest (CoI), DDrift_i and diffusion_i are, respectively, the derivative of the drift (1st Kramers-Moyal term) and the diffusion (2nd Kramers-Moyal term) at the CoI. The first s elements must be monomorphic configurations (only one non-zero entrance, equal to Z) and any positive value for DDrift ad diffusion at those is valid. TransitionMatrix is a sparse array with the transition probabilities between the CoI in the long run. Both ConfigurationsOfInterest and TransitionMatrix can be obtained from TransitionMatrixESML (see ?TransitionMatrixESML for details).";


(* :Code Section: *)
Begin["`Private`"];

Options[TransitionMatrixESML]={"Verbose"->False};
TransitionMatrixESML[s_,T_,Z_,OptionsPattern[]]:=
If[!ValueQ[T[2,1,0]],
Message[TransitionMatrixESML::transitiondefinition,T];
{{{}},{{}}},

With[{prints=OptionValue["Verbose"]},
Module[{CoIFull={{0}},innerCoI={{0.}},\[Rho]Temp={{0}},\[Rho]={{0}},pairId=0,zeros={{0.,0.,0}},TMat,speed={0.},innerCoIcount=0,
CoI=Table[{SparseArray[{i->Z},{s}],1.,1.},{i,1,s}],temp
},


(*****************************************************************************)
If[prints,temp=PrintTemporary["Computing CoI and Transitions"];];
(*****************************************************************************)
{innerCoI,\[Rho]}=Rest[Reap[

Do[pairId++;

{zeros,\[Rho]Temp}=TransitionsPairs[Table[If[i==Z,0.,T[s2,s1,Z-i]],{i,0,Z}],Table[If[i==0,0.,T[s1,s2,i]],{i,0,Z}]];

If[Length[zeros]>2,
Sow[Table[{
SparseArray[{s1->(zeros[[i,1]]-1),s2->(Z-(zeros[[i,1]]-1))},{s}],zeros[[i,3]],zeros[[i,4]]}
,{i,2,Length[zeros]-1}],"CoIinfo"];
,
Sow[{},"CoIinfo"];
];

Do[

Sow[{Which[array[[1,1]]==1,s2,array[[1,1]]==Length[\[Rho]Temp],s1,True,s+innerCoIcount+array[[1,1]]-1],Which[array[[1,2]]==1,s2,array[[1,2]]==Length[\[Rho]Temp],s1,True,s+innerCoIcount+array[[1,2]]-1]}->array[[2]],"traMat"];

,{array,ArrayRules[\[Rho]Temp][[1;;-2]]}];

innerCoIcount+=Length[zeros]-2;

,{s1,1,s-1},{s2,s1+1,s}];
]][[1]];

If[prints,NotebookDelete[temp];];
(*****************************************************************************)

innerCoI=Flatten[innerCoI,1];

If[prints,temp=PrintTemporary["Flatten do CoI"];];
CoI=Flatten[{CoI,innerCoI},1];
If[prints,NotebookDelete[temp];];
(*****************************************************************************)
If[prints,temp=PrintTemporary["Matrix construction"];];
(*****************************************************************************)
TMat=SparseArray[\[Rho],{1,1}(s+innerCoIcount)];

speed=Total[Transpose[TMat]];
TMat/=Max[speed];
speed/=Max[speed]; 



If[prints,NotebookDelete[temp];];
(*****************************************************************************)


Do[TMat[[i,i]]=1-speed[[i]];
,{i,1,Length@TMat}];


{CoI,TMat}

]
]
];


Options[StationaryDistributionESML]={"Verbose"->False,Method->{"Arnoldi","Shift"->1.00001,"Tolerance"-> 0,"MaxIterations"->10^4}};
StationaryDistributionESML[s_,Z_,CoI_,TMat_,opts:OptionsPattern[]]:=
With[{prints=OptionValue["Verbose"]},
Module[{ vec={0.},speed={0.},\[Sigma]2=0.,temp},

If[Length[TMat]>2,
vec=Eigenvectors[Transpose@TMat,1,Evaluate[FilterRules[{opts}, Options[Eigenvectors]]],Method->{"Arnoldi","Shift"->1.00001,"Tolerance"-> 0,"MaxIterations"->10^4}];
,
vec=Eigenvectors[Transpose@TMat,1,Evaluate[FilterRules[{opts}, Options[Eigenvectors]]]];
];

(*****************************************************************************)
If[prints,temp=PrintTemporary["Renormalization"];];
(*****************************************************************************)
Which[
Length[vec]==0,
Print["Could not compute Stationary Distribution."];
vec=ConstantArray[0,Length@CoI];
,
Length[vec]>1,
Print["Stationary Distribution was computed but might contain errors."];
vec=vec[[-1]];
,
True,
vec=vec[[1]];

Do[
If[CoI[[i,2]]<0,
\[Sigma]2=(CoI[[i,3]]//N)/-CoI[[i,2]];
vec[[i]]*=Z Sqrt[\[Pi] \[Sigma]2/2. ](Erf[(Norm[CoI[[i,1]]]-1.)/(Z Sqrt[2. \[Sigma]2])]+Erf[(Z-Norm[CoI[[i,1]]]-1.)/(Z Sqrt[2. \[Sigma]2])])];
,{i,s+1,Length@vec}];

vec=vec/Total[vec];

];

If[prints,NotebookDelete[temp];];


{CoI[[All,1]],vec}


]
];



(* :Error Messages: *)

TransitionMatrixESML::transitiondefinition="The function `1` is undefined. Please guarantee that `1`[s1,s2,i] has definition.";
TransitionMatrixESML::wrongpath="The path `1` is not a valid directory. Please provide a valid directory (e.g., NotebookDirectory[]). Exporting to `2`.";
StationaryDistributionESML::transitiondefinition="The function `1` is undefined. Please guarantee that `1`[s1,s2,i,`2`] has definition.";


(* :Code Section: *)





(* :Code Section: *)
DiscreteZero=Compile[{{vecXY,_Real,2}},
Module[{
dx=Internal`Bag[],
x=Internal`Bag[],
index=Internal`Bag[],
Z=Length[vecXY]-1,
xtemp=0.,
dxtemp=0.,
thereAreZeros=False
},

Do[
If[vecXY[[i,2]]>= 0&&vecXY[[i+1,2]]< 0,
If[vecXY[[i,2]]== 0&&i>1,
Internal`StuffBag[dx,(vecXY[[i+1,2]]-vecXY[[i-1,2]])/(vecXY[[i+1,1]]-vecXY[[i-1,1]])];
Internal`StuffBag[x,vecXY[[i,1]]];
Internal`StuffBag[index,i];
If[!thereAreZeros,thereAreZeros=True];
,
dxtemp=(vecXY[[i+1,2]]-vecXY[[i,2]])/(vecXY[[i+1,1]]-vecXY[[i,1]]);
xtemp=vecXY[[i,1]]-vecXY[[i,2]]/dxtemp; 

If[Abs[xtemp-vecXY[[i,1]]]<Abs[xtemp-vecXY[[i+1,1]]]&&i>1,
Internal`StuffBag[dx,dxtemp];
Internal`StuffBag[x,xtemp];
Internal`StuffBag[index,i+0.];
If[!thereAreZeros,thereAreZeros=True;];
,
If[i<Z,
Internal`StuffBag[dx,dxtemp];
Internal`StuffBag[x,xtemp];
Internal`StuffBag[index,i+1.];
If[!thereAreZeros,thereAreZeros=True;];
];
];
];
];
,{i,1,Z}];



If[thereAreZeros,
Transpose[{Internal`BagPart[dx,All,List],Internal`BagPart[x,All,List],Internal`BagPart[index,All,List]}],
{{}}
]
]
,
{{Internal`BagLength[_],_Integer}},
CompilationTarget->"C"];

TransitionsPairs[Tp_List,Tm_List]:=
If[Length[Tp]==Length[Tm],
With[{Z=Length[Tp]-1},
Module[{zeros=DiscreteZero[Transpose@{Range[0,Z]/Z,Tp-Tm}], CoIindex={1,Z+1},Pp=1.,Pm=1.},

If[zeros!={{}},

(* consider internal transitions *)
CoIindex=Flatten@Insert[CoIindex,IntegerPart[zeros[[All,-1]]],2] ;
{Transpose@{
CoIindex, (* Point index, from 1 to Z+1 *)
Flatten[{0.,zeros[[All,2]],1.}],(*Point position from 0 to 1*)
Flatten[{1.,zeros[[All,1]],1.}],(*Drift derivative*)
Flatten[{1.,Table[(Tm[[CoIindex[[i]]]]+Tp[[CoIindex[[i]]]])/(2Z),{i,2,Length[CoIindex]-1}],1.}](*Diffusion*)
},

(* calculate transitions between COI *)
SparseArray[Flatten[Table[

(*Pp=1+\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(j = CoI[\([\)\(i\)\(]\)] + 1\), \(CoI[\([\)\(i + 1\)\(]\)] - 1\)]\(
\*UnderoverscriptBox[\(\[Product]\), \(k = CoI[\([\)\(i\)\(]\)] + 1\), \(j\)]
\*FractionBox[\(Tm[\([\)\(k\)\(]\)]\), \(Tp[\([\)\(k\)\(]\)]\)]\)\);*)
Pp=1.;Do[Pp=1.+Pp Tm[[k]]/Tp[[k]];,{k,CoIindex[[i+1]]-1,CoIindex[[i]]+1,-1}];
(*Pm=1+\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(j = CoI[\([\)\(i\)\(]\)] + 1\), \(CoI[\([\)\(i + 1\)\(]\)] - 1\)]\(
\*UnderoverscriptBox[\(\[Product]\), \(k = CoI[\([\)\(i\)\(]\)] + 1\), \(j\)]
\*FractionBox[\(Tp[\([\)\(k\)\(]\)]\), \(Tm[\([\)\(k\)\(]\)]\)]\)\)*)
Pm=1.;Do[Pm=1.+Pm Tp[[k]]/Tm[[k]];,{k,CoIindex[[i]]+1,CoIindex[[i+1]]-1}];
(*CoI[[i]] to CoI[[i+1]]  and   CoI[[i+1]] to CoI[[i]]*)

{{i,i+1}-> Tp[[CoIindex[[i]]]] /Pp,{i+1,i}->Tm[[CoIindex[[i+1]]]]/Pm}
,{i,1,Length[CoIindex]-1}]],{Length[CoIindex],Length[CoIindex]}]
},

(* consider transitions between monomorphic states *)
{Transpose@{
CoIindex, (* Point index, from 1 to Z+1 *)
{0.,1.},(*Point position from 0 to 1*)
{1.,1.},(*Drift derivative*)
{1.,1.}(*Diffusion*)
},
SparseArray[

(*Pp=1+\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(j = CoI[\([\)\(i\)\(]\)] + 1\), \(CoI[\([\)\(i + 1\)\(]\)] - 1\)]\(
\*UnderoverscriptBox[\(\[Product]\), \(k = CoI[\([\)\(i\)\(]\)] + 1\), \(j\)]
\*FractionBox[\(Tm[\([\)\(k\)\(]\)]\), \(Tp[\([\)\(k\)\(]\)]\)]\)\);*)
Pp=1.;Do[Pp=1.+Pp Tm[[k]]/Tp[[k]];,{k,2,Z}];
(*Pm=1+\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(j = CoI[\([\)\(i\)\(]\)] + 1\), \(CoI[\([\)\(i + 1\)\(]\)] - 1\)]\(
\*UnderoverscriptBox[\(\[Product]\), \(k = CoI[\([\)\(i\)\(]\)] + 1\), \(j\)]
\*FractionBox[\(Tp[\([\)\(k\)\(]\)]\), \(Tm[\([\)\(k\)\(]\)]\)]\)\)*)
Pm=1.;Do[Pm=1.+Pm Tp[[k]]/Tm[[k]];,{k,Z,2,-1}];
(*CoI[[i]] to CoI[[i+1]]  and   CoI[[i+1]] to CoI[[i]]*)

{{1,2}-> Tp[[1]] /Pp,{2,1}->Tm[[Z+1]]/Pm}
,{Length[CoIindex],Length[CoIindex]}]
}
]

]]
,Message[TransitionsPairs::difdimensions,Dimensions[Tp],Dimensions[Tm]]];



TransitionsPairsSML[Tp_List,Tm_List]:=
If[Length[Tp]==Length[Tm],
With[{Z=Length[Tp]-1,CoIindex={1,Length[Tp]}},
Module[{Pp=1.,Pm=1.},

(* consider transitions between monomorphic states *)
{Transpose@{
CoIindex, (* Point index, from 1 to Z+1 *)
{0.,1.},(*Point position from 0 to 1*)
{1.,1.},(*Drift derivative*)
{1.,1.}(*Diffusion*)
},
SparseArray[

(*Pp=1+\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(j = CoI[\([\)\(i\)\(]\)] + 1\), \(CoI[\([\)\(i + 1\)\(]\)] - 1\)]\(
\*UnderoverscriptBox[\(\[Product]\), \(k = CoI[\([\)\(i\)\(]\)] + 1\), \(j\)]
\*FractionBox[\(Tm[\([\)\(k\)\(]\)]\), \(Tp[\([\)\(k\)\(]\)]\)]\)\);*)
Do[Pp=1.+Pp Tm[[k]]/Tp[[k]];,{k,Z,2,-1}];
(*Pm=1+\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(j = CoI[\([\)\(i\)\(]\)] + 1\), \(CoI[\([\)\(i + 1\)\(]\)] - 1\)]\(
\*UnderoverscriptBox[\(\[Product]\), \(k = CoI[\([\)\(i\)\(]\)] + 1\), \(j\)]
\*FractionBox[\(Tp[\([\)\(k\)\(]\)]\), \(Tm[\([\)\(k\)\(]\)]\)]\)\)*)
Do[Pm=1.+Pm Tp[[k]]/Tm[[k]];,{k,2,Z}];
(*CoI[[i]] to CoI[[i+1]]  and   CoI[[i+1]] to CoI[[i]]*)

{{1,2}->Tp[[1]] /Pp,{2,1}->Tm[[Z+1]]/Pm}
,{2,2}]
}

]]
,Message[TransitionsPairsSML::difdimensions,Dimensions[Tp],Dimensions[Tm]]];

Options[TransitionMatrixSML]={"Verbose"->False};
TransitionMatrixSML[s_,T_,Z_,OptionsPattern[]]:=
If[!ValueQ[T[2,1,0]],
Message[TransitionMatrixESML::transitiondefinition,T];
{{{}},{{}}},

With[{prints=OptionValue["Verbose"]},
Module[{CoIFull={{0}},\[Rho]Temp={{0}},\[Rho]={{0}},zeros={{0.,0.,0}},TMat,speed={0.},
CoI=Table[{SparseArray[{i->Z},{s}],1.,1.},{i,1,s}],temp
},


(*****************************************************************************)
If[prints,temp=PrintTemporary["Computing CoI and Transitions"];];
(*****************************************************************************)
\[Rho]=Flatten@
Table[
{zeros,\[Rho]Temp}=TransitionsPairsSML[Table[If[i==Z,0.,T[s2,s1,Z-i]],{i,0,Z}],Table[If[i==0,0.,T[s1,s2,i]],{i,0,Z}]];
Table[
{If[array[[1,1]]==1,s2,s1],If[array[[1,2]]==1,s2,s1]}->array[[2]]
,{array,Most[ArrayRules[\[Rho]Temp]]}]
,{s1,1,s-1},{s2,s1+1,s}];

If[prints,NotebookDelete[temp];];
(*****************************************************************************)


(*****************************************************************************)
If[prints,temp=PrintTemporary["Matrix construction"];];
(*****************************************************************************)
TMat=SparseArray[\[Rho],{s,s}];

speed=Total[Transpose[TMat]];
(*TMat/=Max[speed];
speed/=Max[speed];*)


If[prints,NotebookDelete[temp];];
(*****************************************************************************)


Do[TMat[[i,i]]=1-speed[[i]];
,{i,1,Length@TMat}];

{CoI,TMat}

]
]
];



End[];


(* :Code Section (Call Protect): *)


Protect[TransitionMatrixESML];
Protect[StationaryDistributionESML];
Protect[DiscreteZero];
Protect[TransitionsPairs];
Protect[TransitionsPairsSML];
Protect[TransitionMatrixSML];





EndPackage[];
