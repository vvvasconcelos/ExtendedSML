(* ::Package:: *)

DiscreteZero[vecxy_List]:=
Module[{vecsi,vecsi2,vvecsi,vecy,transi,x,dx,ToConsider,CloseIndex,dxtemp=0.,xtemp=0.},
  vecy=vecxy[[All,2]];(*list of y values*)
  vecsi=Transpose@{Sign[vecy],Range[Length[vecy]]};(*List of signal (-1,0,1) for each index*)
  vecsi2=Select[vecsi,First[#]!=0&];(*Filter out zeros: List of signal (-1,1) for each index*)
  vvecsi=Split[vecsi2,First[#1]==First[#2]&];(*Group into positive values of y and negative values*)
  transi=Transpose@{Most[Max[#[[All,2]]]&/@vvecsi],Rest[Min[#[[All,2]]]&/@vvecsi]};(*index of transitions from positive to negative values (or vice versa)*)
  ToConsider=Table[If[vecsi[[transi[[i,1]],1]]>0,True,False],{i,1,Length[transi]}];(*Check if it is a stable point of not: List of True and False for each transition*)

  {dx,x,CloseIndex}=
    Transpose@
      Table[If[ToConsider[[i]], 
        {dxtemp=(vecxy[[transi[[i,2]],2]]-vecxy[[transi[[i,1]],2]])/(vecxy[[transi[[i,2]],1]]-vecxy[[transi[[i,1]],1]]),
          xtemp=vecxy[[transi[[i,1]],1]]-vecxy[[transi[[i,1]],2]]dxtemp,
          Ordering[Abs[vecxy[[All,1]]-xtemp],1][[1]]}
          ,{0.,0.,0}]
         ,{i,1,Length[transi]}
      ];

  Select[Select[Transpose@{x,dx,CloseIndex,ToConsider},#[[-1]]&][[All,1;;-2]],#[[-1]]!=1&&#[[-1]]!=Length[vecy]&]
];


DiscreteZero2[vecxy_List]:=With[{init=10},
Module[{x=Table[0.,{i,1,init}],dx=Table[0.,{i,1,init}],CloseIndex=Table[0,{i,1,init}],vsize=init,nzeros=0,dxtemp=0.,xtemp=0.},


Do[
If[vecxy[[i,2]]> 0&&vecxy[[i+1,2]]< 0,
If[nzeros+1>vsize,
AppendTo[x,Table[0.,{i,1,init}]];
AppendTo[dx,Table[0.,{i,1,init}]];
AppendTo[CloseIndex,Table[0,{i,1,init}]];
vsize+=init;
];
nzeros++;

dx[[nzeros]]=(vecxy[[i+1,2]]-vecxy[[i,2]])/(vecxy[[i+1,1]]-vecxy[[i,1]]);
x[[nzeros]]=vecxy[[i,1]]-vecxy[[i,2]]dx[[nzeros]];
CloseIndex[[nzeros]]=If[Abs[xtemp-vecxy[[i,1]]]<Abs[xtemp-vecxy[[i+1,1]]],
i,i+1
];

];
,{i,1,Length[vecxy]-1}];

Transpose@{x[[1;;nzeros]],dx[[1;;nzeros]],CloseIndex[[1;;nzeros]]}
]];
