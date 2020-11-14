(* ::Package:: *)

BeginPackage["ReadAMBiTOutput`"];


ReadSingleTransition::usage = "Parses AMBiT output text file and creates list of Single Shell Orbital E1 Transitions"

ReadOrbitals::usage = "Parses AMBiT output text file and creates list of Single Shell Orbitals and their energiers in a.u."

ReadCI::usage = "Parses AMBiT output text file and creates lists for all CI levels. Returns lists in form {PI,{J,Pi,N,{Ind,E,{Config,Percent}}}}. Pi=even->1, Pi=odd->2."

DecomposeTransitions::usage = "DecomposeTransitions[STR] splits a line of the transition list into its constituent parts. Iterate of this function in a table. Returns {Jinitial,Pinitial,CIinitial,Jfinal,Pfinal,CIfinal,S}."

ReadTransitions::usage = "Parses AMBiT output text file and creates lists for E1 CI transition. Returns lists in form {Jinitial,Pinitial,CIinitial,Jfinal,Pfinal,CIfinal,S}"

ReadLevels::usage = "Parses AMBiT output text file and creates lists for Average Configuration Levels. Returns lists for {Configurations,Energy,Num}"*)


Begin["`Private`"];


ReadSingleTransition[FILE_]:=Module[{SingleShellPos,ShellPosStart,ShellPosEnd,i,COND,STR,Single,SingleSplit},
{
	SingleShellPos=Position[FILE,"One-body transition reduced matrix elements (a.u.):"];
	ShellPosStart=SingleShellPos[[1]][[1]];
	i=ShellPosStart;
	COND=True;
	While[COND==True,{
		STR=FILE[[i+1]];
		If[STR=={},COND=False,i=i+1];
	}];
	ShellPosEnd=i;
	Single=FILE[[ShellPosStart+1;;ShellPosEnd]];
	L=Length[Single];
	SingleSplit={};
	For[i=1,i<=L,i=i+1,{
		STR=StringSplit[Single[[i]]][[1]];
		AppendTo[SingleSplit,STR];
	}];
};
Return[SingleSplit];
];

ReadOrbitals[FILE_]:=Module[{OrbitalPosStart,COND,i,ExcitedPos,STR,OrbitalPosEnd,Core,Excited,HF,L,OrbitalEnergy,Orbitals,Energy},
{
	OrbitalPosStart=Position[FILE,"Core orbitals:"][[1]][[1]];
	COND=True;
	i=OrbitalPosStart;
	ExcitedPos=0;
	While[COND==True,
	{
		STR=FILE[[i+1]];
		If[STR=={},COND=False,
		{
			If[StringMatchQ[STR[[1]],"Excited orbitals:"]==True,ExcitedPos=i+1];
			i=i+1;
		}];
	}];
	OrbitalPosEnd=i;
	Core=FILE[[OrbitalPosStart+1;;ExcitedPos-1]];
	Excited=FILE[[ExcitedPos+1;;OrbitalPosEnd]];
	HF=Join[Core,Excited];
	HF=Union[HF];
	L=Length[HF];
	
	OrbitalEnergy=Table[StringSplit[HF[[i]]][[1]],{i,1,L}];
	
	Orbitals=OrbitalEnergy\[Transpose][[1]];
	Energy=ToExpression[OrbitalEnergy\[Transpose][[4]]];
};
Return[{Orbitals,Energy}];
];




ReadCI[FILE_]:=Module[{PosJ,LJ,CIsol,STR,STRSp,x,J,P,Num,COND,i,j,CIJ,Configuration,Percentage,CHR,Config,Perc,CI1,CHR1,a,Ind,Energy},
{
	PosJ=Position[FILE,s_String/;StringMatchQ[s,"Solutions for J*"]];
	LJ=Length[PosJ];
	CIsol=ConstantArray[{},2];
	For[j=1,j<=LJ,j=j+1,
	{
		STR=FILE[[PosJ[[j]][[1]]]];
		STRSp=StringSplit[STR];
		x=Position[STRSp,"="];
		J=ToExpression[STRSp[[x[[1]][[1]]]][[x[[1]][[2]]+1]]];
		P=STRSp[[x[[2]][[1]]]][[x[[2]][[2]]+1]];
		If[P=="even",P=1,P=2];
		Num=STRSp[[x[[3]][[1]]]][[x[[3]][[2]]+1]];
		Num=ToExpression[StringTake[Num,{1,-3}]];

		COND=True;
		i=PosJ[[j]][[1]];
		CIJ={};
		Configuration={};
		Percentage={};
		While[COND==True, 
		{
			i=i+1;
			STR=FILE[[i]];
			If[STR!={},
			{
				If[StringFreeQ[STR,"Number"][[1]]==True 
				&& StringFreeQ[STR,"E1 transition strengths (S):"][[1]]==True 
				&& StringFreeQ[STR,"Solutions"][[1]]==True
				&& StringFreeQ[STR,"Finished"][[1]]==True,
				{
					CHR=StringSplit[STR,"  "];
					If[Length[CHR[[1]]]==2,
					{
						Config=CHR[[1]][[1]];
						Configuration=Append[Configuration,Config];
						Perc=Internal`StringToDouble[CHR[[1]][[2]]]*0.01;
						Percentage=Append[Percentage,Perc];
					},{
						CHR1=CHR[[1]][[1]];
						a=StringPosition[CHR1,":"];
						Ind=StringTake[CHR1,{1,a[[1]][[1]]-1}];
						Energy=ToExpression[StringTake[CHR1,{a[[1]][[1]]+1,-1}]];
					}];
				}, {
				COND=False
				}];
			},{
				CI1={Ind,Energy,Configuration,Percentage};
				Configuration={};
				Percentage={};
				If[CI1[[4]]=={},Continue[],CIJ=Append[CIJ,CI1]];
			}];
		}];
		AppendTo[CIsol[[P]],{J,P,Num,CIJ}];
	}];
};
Return[CIsol];
];



Parity=<|"e"->1,"o"->2|>;
DecomposeTransitions[STR_]:=Module[{P1,P2,Jinitial,Jfinal,Pinitial,Pfinal,CIinitial,CIfinal,S,CHR,s,Initial,a,Final},
{
	CHR=StringSplit[STR];
	s=CHR[[1]][[-1]];
	S=Internal`StringToDouble[s];

	Initial=CHR[[1]][[1]];
	a=StringPosition[Initial,":"];
	Jinitial=Internal`StringToDouble[StringTake[Initial,{1,a[[1]][[1]]-2}]]/2//N;
	P1=StringTake[Initial,{a[[1]][[1]]-1}];
	(*If[P1=="e",Pinitial=1,Pinitial=2];*)
	Pinitial=Parity[[P1]];
	CIinitial=IntegerPart[Internal`StringToDouble[StringTake[Initial,{a[[1]][[1]]+1,-1}]]];

	Final=CHR[[1]][[3]];
	a=StringPosition[Final,":"];
	Jfinal=Internal`StringToDouble[StringTake[Final,{1,a[[1]][[1]]-2}]]/2//N;
	P2=StringTake[Final,{a[[1]][[1]]-1}];
	(*If[P2=="e",Pfinal=1,Pfinal=2];*)
	Pfinal=Parity[[P2]];
	CIfinal=IntegerPart[Internal`StringToDouble[StringTake[Final,{a[[1]][[1]]+1,-1}]]];
};
Return[{Jinitial,Pinitial,CIinitial,Jfinal,Pfinal,CIfinal,S}]
]

ReadTransitions[FILE_]:=Module[{TransitionPos,TPosStart,TPosEnd,Transitions,L,T},
{
	TransitionPos=Position[FILE,"E1 transition strengths (S):"];
	If[Length[TransitionPos]==3,
	{
		TransitionPos=TransitionPos[[2;;3]];
	}];
	TPosStart=TransitionPos[[1]][[1]];
	TPosEnd=TransitionPos[[2]][[1]]-2;

	Transitions=FILE[[TPosStart+1;;TPosEnd]];

	L=Length[Transitions];
	(*Tic=AbsoluteTime[];*)
	T=Table[DecomposeTransitions[Transitions[[i]]],{i,1,L}];
(*	Toc=AbsoluteTime[];
	Print[L," -> ",Toc-Tic];
	Print[1," -> ",(Toc-Tic)/L];*)
};
Return[T\[Transpose]];
];


ReadLevels[FILE_]:=Module[
{LevelsPosStart,COND,STR,LevelsPosEnd,L,Levels,Levels1,Shell,Energy,Num,ord,i},
{
	LevelsPosStart=Position[FILE,"Configurations; E; N:"][[1]][[1]];
	i=LevelsPosStart;
	While[FILE[[i]]!={},
	{
		i=i+1;
	}];
	LevelsPosEnd=i-1;
	Levels=FILE[[LevelsPosStart+1;;LevelsPosEnd]];
	L=Length[Levels];
	Levels1=Table[StringSplit[Levels[[i]],";"][[1]],{i,1,L}];

	Shell=Levels1\[Transpose][[1]];
	Energy=ToExpression[Levels1\[Transpose][[2]]];
	Num=ToExpression[Levels1\[Transpose][[3]]];
	ord=Ordering[Energy];
	Shell=Shell[[ord]];
	Energy=Energy[[ord]];
	Num=Num[[ord]];

};
Return[{Shell,Energy,Num}];
];


End[];
EndPackage[];
