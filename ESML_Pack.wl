(* ::Package:: *)

BeginPackage["ESMLPack`"];


StatDESML::usage="StatDESML[s,Z,T,{par1,par2...}] is this and that and gives that and that";


Begin["Private`"];


Needs["CCompilerDriver`"];
If[Length[CCompilers[]]==0,
$CCompiler={"Name"->"Intel Compiler","Compiler"->CCompilerDriver`IntelCompiler`IntelCompiler,"CompilerInstallation"->"D:\\Program Files (x86)\\Intel\\Composer XE 2015","CompilerName"->Automatic};
];
DiscreteZero=Compile[{{vecXY,_Real,2}},
Module[{
dx=Internal`Bag[Most[{0.}]],
x=Internal`Bag[Most[{0.}]],
index=Internal`Bag[Most[{0}]],
Z=Length[vecXY]-1,
xtemp=0.,
dxtemp=0.,
thereAreZeros=False
},

Do[
If[vecXY[[i,2]]>= 0&&vecXY[[i+1,2]]< 0,
If[vecXY[[i,2]]== 0&&i>1,
Internal`StuffBag[dx,(vecXY[[i+1,2]]-If[i>1,vecXY[[i-1,2]],vecXY[[i,2]]])/(vecXY[[i+1,1]]-If[i>1,vecXY[[i-1,1]],vecXY[[i,1]]])];
Internal`StuffBag[x,vecXY[[i,1]]];
Internal`StuffBag[index,i];
If[!thereAreZeros,thereAreZeros=True;];
,
dxtemp=(vecXY[[i+1,2]]-vecXY[[i,2]])/(vecXY[[i+1,1]]-vecXY[[i,1]]);
xtemp=vecXY[[i,1]]-vecXY[[i,2]]/((vecXY[[i+1,2]]-vecXY[[i,2]])/(vecXY[[i+1,1]]-vecXY[[i,1]])); 

If[Abs[xtemp-vecXY[[i,1]]]<Abs[xtemp-vecXY[[i+1,1]]]&&i>1,
Internal`StuffBag[dx,dxtemp];
Internal`StuffBag[x,xtemp];
Internal`StuffBag[index,i];
If[!thereAreZeros,thereAreZeros=True;];
,
If[i>1&&i<Z,
Internal`StuffBag[dx,dxtemp];
Internal`StuffBag[x,xtemp];
Internal`StuffBag[index,i+1];
If[!thereAreZeros,thereAreZeros=True;];
];
];];
];
,{i,1,Z}];



If[thereAreZeros,
Transpose[{Internal`BagPart[dx,All,List],Internal`BagPart[x,All,List],Internal`BagPart[index,All,List]}],
{{}}
]
]
,CompilationTarget->"C"];


TransitionsPairs::difdimensions="The lenght of the first list, `1`, is equal to that of the second, `2`.";
TransitionsPairs[Tp_List,Tm_List]:=If[Length[Tp]==Length[Tm],With[{Z=Length[Tp]-1},
Module[{zeros=DiscreteZero[Transpose@{Range[0,Z]/Z,Tp-Tm}], CoIindex={1,Z+1},Pp=1.,Pm=1.},

If[zeros!={{}},
CoIindex=Flatten@Insert[CoIindex,IntegerPart[zeros[[All,-1]]],2] ;


{Transpose@{
CoIindex, (* Point index, from 1 to Z+1 *)
Flatten[{0.,zeros[[All,2]],1.}],(*Point position from 0 to 1*)
Flatten[{1.,zeros[[All,1]],1.}],(*Drift derivative*)
Flatten[{1.,Table[(Tm[[CoIindex[[i]]]]+Tp[[CoIindex[[i]]]])/(2Z),{i,2,Length[CoIindex]-1}],1.}](*Diffusion*)
},
SparseArray[Flatten[Table[

(*Pp=1+\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(j = CoI[\([i]\)] + 1\), \(CoI[\([i + 1]\)] - 1\)]\(
\*UnderoverscriptBox[\(\[Product]\), \(k = CoI[\([i]\)] + 1\), \(j\)]
\*FractionBox[\(Tm[\([k]\)]\), \(Tp[\([k]\)]\)]\)\);*)
Pp=1.;Do[Pp=1.+Pp Tm[[k]]/Tp[[k]];,{k,CoIindex[[i]]+1,CoIindex[[i+1]]-1}];
(*Pm=1+\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(j = CoI[\([i]\)] + 1\), \(CoI[\([i + 1]\)] - 1\)]\(
\*UnderoverscriptBox[\(\[Product]\), \(k = CoI[\([i]\)] + 1\), \(j\)]
\*FractionBox[\(Tp[\([k]\)]\), \(Tm[\([k]\)]\)]\)\)*)
Pm=1.;Do[Pm=1.+Pm Tp[[k]]/Tm[[k]];,{k,CoIindex[[i+1]]-1,CoIindex[[i]]+1,-1}];
(*CoI[[i]] to CoI[[i+1]]  and   CoI[[i+1]] to CoI[[i]]*)

{{i,i+1}-> Tp[[CoIindex[[i]]]] /Pp,{i+1,i}->Tm[[CoIindex[[i+1]]]]/Pm}
,{i,1,Length[CoIindex]-1}]],{Length[CoIindex],Length[CoIindex]}]
},


{Transpose@{
CoIindex, (* Point index, from 1 to Z+1 *)
{0.,1.},(*Point position from 0 to 1*)
{1.,1.},(*Drift derivative*)
{1.,1.}(*Diffusion*)
},
SparseArray[

(*Pp=1+\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(j = CoI[\([i]\)] + 1\), \(CoI[\([i + 1]\)] - 1\)]\(
\*UnderoverscriptBox[\(\[Product]\), \(k = CoI[\([i]\)] + 1\), \(j\)]
\*FractionBox[\(Tm[\([k]\)]\), \(Tp[\([k]\)]\)]\)\);*)
Pp=1.;Do[Pp=1.+Pp Tm[[k]]/Tp[[k]];,{k,2,Z}];
(*Pm=1+\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(j = CoI[\([i]\)] + 1\), \(CoI[\([i + 1]\)] - 1\)]\(
\*UnderoverscriptBox[\(\[Product]\), \(k = CoI[\([i]\)] + 1\), \(j\)]
\*FractionBox[\(Tp[\([k]\)]\), \(Tm[\([k]\)]\)]\)\)*)
Pm=1.;Do[Pm=1.+Pm Tp[[k]]/Tm[[k]];,{k,Z,2,-1}];
(*CoI[[i]] to CoI[[i+1]]  and   CoI[[i+1]] to CoI[[i]]*)

{{1,2}-> Tp[[1]] /Pp,{2,1}->Tm[[Z+1]]/Pm}
,{Length[CoIindex],Length[CoIindex]}]
}


]

]]
,Message[TransitionsPairs::difdimensions,Dimensions[Tp],Dimensions[Tm]]];


StatDESML::transitiondefinition="The function `1` is undefined. Please guarantee that `1`[s1,s2,i,`2`,`3`] has definition.";
StatDESML[s_,T_,Z_,par_List]:=If[ValueQ[T[2,1,0,Z,par]],
With[{prints=False},
Module[{CoIFull={{0}},innerCoI={{0.}},\[Rho]Temp={{0}},\[Rho]={{0}},pairId=0,zeros={{0.,0.,0}},TMat, vec={0.},speed={0.},temp
,innerCoIcount=0,
CoI=Table[{SparseArray[{i->Z},{s}],1.,1.},{i,1,s}]
},







(*****************************************************************************)
If[prints,temp=PrintTemporary["Computing CoI and Transitions"];];
(*****************************************************************************)
{innerCoI,\[Rho]}=Rest[Reap[

Do[pairId++;

{zeros,\[Rho]Temp}=TransitionsPairs[Table[T[s2,s1,Z-i,Z,par],{i,0,Z}],Table[T[s1,s2,i,Z,par],{i,0,Z}]];

If[Length[zeros]>2,
Sow[Table[{
SparseArray[{s1->(zeros[[i,1]]-1),s2->(Z-(zeros[[i,1]]-1))},{s}],zeros[[i,3]],zeros[[i,4]]}
,{i,2,Length[zeros]-1}],CoIinfo];
,
Sow[{},CoIinfo];
];

Do[

Sow[{Which[array[[1,1]]==1,s2,array[[1,1]]==Length[\[Rho]Temp],s1,True,s+innerCoIcount+array[[1,1]]-1],Which[array[[1,2]]==1,s2,array[[1,2]]==Length[\[Rho]Temp],s1,True,s+innerCoIcount+array[[1,2]]-1]}->array[[2]],traMat];

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

speed=Total[Transpose@TMat];
TMat/=Max[speed];
speed/=Max[speed]; 


Do[
TMat[[i,i]]=1-speed[[i]];
,{i,1,Length@TMat}];

If[prints,NotebookDelete[temp];];
(*****************************************************************************)

(*****************************************************************************)
If[prints,temp=PrintTemporary["Stat Dist"];];
(*****************************************************************************)
(*vec=NullSpace[Transpose@TMat];  Cannot use for large matrix due to memory constrains. It is more accurate.*)
vec=Eigenvectors[Transpose@TMat,1,Method->{"Arnoldi","Shift"->1.00000000001}];

If[prints,NotebookDelete[temp];];
(*****************************************************************************)

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
If[CoI[[i,2]]<0,vec[[i]]*=Z Sqrt[2\[Pi] (CoI[[i,3]]/-CoI[[i,2]]) ]];
,{i,s+1,Length@vec}];

vec=vec/Total[vec];

];

If[prints,NotebookDelete[temp];];


{CoI[[All,1]],vec}


]
],

Message[StatDESML::transitiondefinition,T,Z,par];
{{},{}}
];


End[];


EndPackage[];
