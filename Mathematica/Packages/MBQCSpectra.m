(* ::Package:: *)

BeginPackage["MBQCSpectra`"];


LevelDensity::usage = "LevelDensity[Energy,Num,G,increment,max,min,Cutoff] returns lists {z, Density}. Density is energy level density obtained from convoluting the configurations with Lorenzian of spreading width G and cutoff energy Cutoff. z=Range[min,max,increment] over which Density is calculated."

LevelSpacing::usage = "LevelSpacing[PsiE,Energy,range] returns Spacing. Takes the list of energies from ComplexWaveEnergies PSiE and finds the level spacing by looking in the interval range centered at Energy"

ComplexWaveEnergies::usage = "ComplexWaveEnergies[G,z,Density,increment,EnergyCI,Num] returns {Emin,PsiE}. Calculates the energies PsiE of all individual states based on the Density obtained from LevelDensity"

ComplexCoefficients::usage = "ComplexCoefficients[G,EnergyCI,Num,Cutoff] returns Components. Gives the eigenstate components in the single particle configuration basis"

ComplexCoefficients2::usage = "ComplexCoefficients[G,EnergyCI,EnergyCSF,Num,Cutoff] returns Components. Gives the eigenstate components in the single particle configuration basis when using average energies that have been shifted due to CI."

ShellList::usage = "ShellList[Config] returns an ordered list of all single electron shells including repetitions in Config."

ShellDifferences::usage = "ShellDifferences[LIST1,LIST2] returns {Diff1,Diff2}. Finds and returns the differing elements with repetition between LIST1 and LIST2, LIST obtained from ShellList."

SingleElectronDipoleSymetric::usage = "SingleElectronDipoleSymetric[DiffInitial,DiffFinal,Transitions] returns the line strength S for the relevant single electron transition. DiffInitial and DiffFinal are lists of differences obtained from ShellDifferences. Transitions is list of single shell transitions obtained from ReadSingle Transition in ReadAMBiTOutput Package."

Jsingle::usage = "Jsingle[SHELL] returns the numerical momentum j from the relativisitic single electron orbital symbol in the SHELL (i.e. the Diff list from ShellDifferences)."

NonRelativisticConfiguration::usage = "NonRelativisticConfiguration[Shells,Num] returns {CONFIGS, NumPROJECT, NonRelConfigs, NumCONFIS, RelConfigs}. CONFIGS and NumPROJECT are the unique nonrealtivistic configurations included and their total number of projections. NonRelConfigs gives the nonrelativistic configuration for each relativistic inclusion in AMBiT. NumCONFIGS gives the number of Rel configs in each NonRel config. RelConfigs gives the relativisitc configuration for each configuration. Ordering of orbitals in the configurations is given by RelOrbitals and NonRelOrbitals in initialisation."

ConigurationAverageJ::usage = "ConigurationAverageJ[CONFIGS] returns the average total momentum J for the list of non-relativistic configurations CONFIGS by computing all possible combinations."

DiscreteLevels::usage = "DiscreteLevels[FILE,Temp] returns {Egap,S,EINITIAL,EFINAL,Intensity,Einstein}. Reads CI and CI Transitions from AMBiT output FILE and calculates list of transition energies Egap, line strengths S, intitial and final energies EINITIAL and EFINAL of the transitions, the transition Intensity at a given Temp, and the Einstein coefficient. Note transitions are given for the absorption spectrum; switch initial and final for emission spectrum."

DiscreteSpectrum::usage = "DiscreteSpectrum[Egap,S,w,G] returns the spectrum from DiscreteLevels convoluted with Lorentzian, of width G over range of photon energies w."

SingleTransitionData::usage = "SingleTransitionData[FILE] take the AMBiT output file and returns list of single electron transition properties T calculated from ComposeSingle"

ComposeSingle::usage = "ComposeSingle[STR,Orbitals,EnergySingle] returns {Dmel,DeltaOmega,ShellInitial,ShellFinal,jInitial,nInitialPos,jFinal,nFinalPos}. Takes a line of the AMBiT output STR, all single electron orbitals Orbitals and energies EnergySingle from AMBiT, and returns the transitions properties. Dmel is the dipole matrix element, DeltaOmega is the energy difference. Shell is orbital, j is momentum"

E1Transitions::usage = "E1Transitions[i,j,Shells,T,nAverage] returns {E1,exp}. Determines whether a single electron transition can occur between two states. Iterates over all states in a Table {i,j}. Takes list of configurations Shells, SingleTransitionData T, and the average occupation of each configuration nAverage and returns index of T for the transition that occurs between states i and j. Also returns the expectation value exp."

ExpVal::usage = "ExpVal[k,Ck2,E1] extracts the expectation value calculated in E1Transitions table E1 for transition k and weighted by eigenstate coefficient Ck2"

ListIndicies::usage = "ListIndicies[L] returns the single dimension lists {iPsi,jPsi} for the initial and final configuration indicies such that the energy of jPsi>iPsi."

Main::usage = "Main[ind,Cutoff,Renorm,G,T,Light,indicies,Dspacing2,Jin,NumLines,E1,Ltran,Components2] iterates through each set of indicies ind in a table to return the line strength for a transition"


Begin["`Private`"];
<<ReadAMBiTOutput`

(* lists used in multiple functions*)
shell={"s","p","p+","d","d+","f","f+","g","g+"}; 
RelOrbitals={"4s","5s","4p","4p+","5p","5p+","4d","4d+","5d","5d+","4f","4f+","5f","5f+","5g","5g+"};
(* s must ne first; nj, nj+ must be consecutive *) 
NonRelOrbitals={"4s","5s","4p","5p","4d","5d","4f","5f","5g"};
l={0,0,1,1,2,2,3,3,4};
MaxElectrons=2(2l+1);
spin={0.5,-0.5};
eV=27.211;


LevelDensity[Energy_,Num_,G_,increment_,max_,min_,Cutoff_]:=Module[{Renorm,pos,z,Density,g,i,MIN},
{
	MIN=Floor[min]; (*rounds down to nearest integer for consistent position finding later *)
	z=Range[MIN,max,increment]; (* Range of energies *)
	z=Round[z,increment]; (* Round to account for machine precision errors *)
	g=G/2; (* width *)
	Renorm=2/\[Pi]*ArcTan[Cutoff/g];
	Density=Sum[
		1/(\[Pi] g) g^2/(g^2+(z-Energy[[i]])^2) Num[[i]]/Renorm * HeavisideTheta[Cutoff-(z-Energy[[i]])]*HeavisideTheta[Cutoff+(z-Energy[[i]])]
		,{i,1,Length[Energy]},Method->"Procedural"];
};
Return[{z,Density}]
];

LevelSpacing[PsiE_,Energy_,range_]:=Module[{Near,Spacing},
{
	Near=Select[PsiE,Energy-range/2<#<Energy+range/2&];
	If[Length[Near]==1,
	{
		Spacing=LevelSpacing[PsiE,Energy,range+1];
	},{
		Spacing=Mean[Differences[Near]];
	}];
};
Return[Spacing];
];


ComplexWaveEnergies[G_,z_,Density_,increment_,EnergyCI_,Num_]:=Module[{k,F,inc,L,NumTotal,PsiBounds,den,LMax,j,i,PsiE,Emin,Components,Ck},
{
	NumTotal=Total[Num];
	PsiBounds=ConstantArray[0,NumTotal+1];
	den=0;
	L=Length[EnergyCI];
	LMax=Length[Density];
	j=1;
	While[Density[[j]]==0,
	{
		j=j+1;
	}];
	PsiBounds[[1]]=z[[j-1]];
	For[i=2,i<=NumTotal+1 ,i=i+1,
	{
		While[den<=1,
		{
			If[j>LMax,
			{
				Print["Error i=",i," last=",PsiBounds[[i-1]]];
				Break[];
			}];
			den=den+Density[[j]]*increment;
			j=j+1;
		}];
		F=Floor[den];
		If[F==1,
		{
			PsiBounds[[i]]=z[[j-1]];
		},{
			If[F!= 0,
			{
				inc=increment/F;
				For[k=1,k<=F,k=k+1,
				{
					PsiBounds[[i]]=k*inc+z[[j-2]];
					i=i+1;
				}];
				i=i-1;
			},{
				PsiBounds[[i]]=PsiBounds[[i-1]]+increment;
			}];
		}];
		den=den-F;
		If[Mod[i,200000]==0,Print[i," of ",NumTotal]];
	}];

	PsiE=(PsiBounds[[1;;-2]]+PsiBounds[[2;;-1]])/2;
	Emin=Min[PsiE];

};
Return[{Emin,PsiE}];
]


ComplexCoefficients[G_,EnergyCI_,Num_,Cutoff_]:=Module[{L,j,i,Components,Ck,Renorm},
{
	L=Length[EnergyCI];
	Renorm=2/\[Pi]*ArcTan[Cutoff/(G/2)];
	Components=ConstantArray[0,{L,L}];
	For[i=1,i<=L,i=i+1,
	{
		Ck=ConstantArray[0,L];
		For[j=1,j<=L,j=j+1,
		{
			Ck[[j]]=Num[[j]]*(G/2)/((EnergyCI[[i]]-EnergyCI[[j]])^2+(G/2)^2)*Renorm*
			HeavisideTheta[Cutoff-(EnergyCI[[i]]-EnergyCI[[j]])]*HeavisideTheta[Cutoff+(EnergyCI[[i]]-EnergyCI[[j]])];
		}];
		Ck=Normalize[Ck];
		Components[[i]]=Ck;
	}]
};
Return[Components];
]

ComplexCoefficients2[G_,EnergyCI_,EnergyCSF_,Num_,Cutoff_]:=Module[{L,j,i,Components,Ck,Renorm},
{
	L=Length[EnergyCI];
	Renorm=2/\[Pi]*ArcTan[Cutoff/(G/2)];
	Components=ConstantArray[0,{L,L}];
	For[i=1,i<=L,i=i+1,
	{
		Ck=ConstantArray[0,L];
		For[j=1,j<=L,j=j+1,
		{
			Ck[[j]]=Num[[j]]*(G/2)/((EnergyCI[[i]]-EnergyCSF[[j]])^2+(G/2)^2)*Renorm*
			HeavisideTheta[Cutoff-(EnergyCI[[i]]-EnergyCSF[[j]])]*HeavisideTheta[Cutoff+(EnergyCI[[i]]-EnergyCSF[[j]])];
		}];
		Ck=Normalize[Ck];
		Components[[i]]=Ck;
	}]
};
Return[Components];
]


ShellList[Config_]:=Module[{S,L,list,i,str,x,y,M,SH,LI},
{
	S=StringSplit[Config];
	L=S//Length;
	list={};
	For[i=1,i<=L,i++,
	{
		str=S[[i]];
		x=StringPosition[str,shell]; (* shell is Global list of l used in AMBiT *)
		y=x[[-1]][[2]];
		M=ToExpression[StringTake[str,{y+1,-1}]];
		SH=StringTake[str,{1,y}];
		LI=ConstantArray[SH,M];
		list=Join[list,LI];
	}];
	list=list[[Ordering[list]]];
};
Return[list]
];

ShellDifferences[LIST1_,LIST2_]:=Module[{L,i,j,Diff1,Diff2,Compare},
{
	L=Length[LIST1];
	i=1;j=1;
	Diff1={};
	Diff2={};
	While[i<=L&&j<=L,
	{
		Compare=Order[LIST1[[i]],LIST2[[j]]]; (* compares the ordering of the shells *)
		If[Compare==0,
		{
			i=i+1;
			j=j+1;
		},{
			If[Compare==1,
			{
				AppendTo[Diff1,LIST1[[i]]];
				i=i+1;
			},{
				AppendTo[Diff2,LIST2[[j]]];
				j=j+1;
			}];
		}]
	}];
	If[i!=L+1, (* Process for when 1 of the LISTS has reached the end *)
	{
		If[Diff1=={},Diff1=LIST1[[i;;-1]],Join[Diff1,LIST1[[i;;-1]]]];
	}];
	If[j!=L+1,
	{
		If[Diff2=={},Diff2=LIST2[[j;;-1]],Join[Diff2,LIST2[[j;;-1]]]];
	}];
};
Return[{Diff1,Diff2}];
];

SingleElectronDipoleSymetric[DiffInitial_,DiffFinal_,Transitions_]:=Module[{El,DipMel,S},
{
	El=Select[Transitions,#[[1]]==DiffInitial[[1]]&& #[[3]]==DiffFinal[[1]]&];
	If[El=={},
	{
		El=Select[Transitions,#[[3]]==DiffInitial[[1]]&& #[[1]]==DiffFinal[[1]]&];
	}];
	If[El=={},
	{
		S=0;
	},{
		DipMel=El[[1]][[-1]];
		S=Abs[ToExpression[DipMel]]^2
	}];
};
Return[S]
];

Jsingle[SHELL_]:=Module[{jstr,j,x},
{
	x=StringPosition[SHELL[[1]],shell];
	jstr=StringTake[SHELL[[1]],x[[-1]]];
	Which[
		jstr=="s",j=0.5,
		jstr=="p",j=0.5,
		jstr=="p+",j=1.5,
		jstr=="d",j=1.5,
		jstr=="d+",j=2.5,
		jstr=="f",j=2.5,
		jstr=="f+",j=3.5,
		jstr=="g",j=3.5,
		jstr=="g+",j=4.5
	];
};
Return[j]
];

shell={"s","p","p+","d","d+","f","f+","g","g+"}; 
RelOrbitals={"4s","5s","4p","4p+","5p","5p+","4d","4d+","5d","5d+","4f","4f+","5f","5f+","5g","5g+"};
(* s must ne first; nj, nj+ must be consecutive *) 
NonRelOrbitals={"4s","5s","4p","5p","4d","5d","4f","5f","5g"};
l={0,0,1,1,2,2,3,3,4};

NonRelativisticConfiguration[Shells_,Num_]:=Module[{Lno,Lro,RelConfigs,RelOrbits,NonRelConfigs,list,NonRelOrbits,n,n1,n2,i,j,k,CONFIGS,NumPROJECT,NumCONFIGS,pos},
{
	Lro=Length[RelOrbitals];
	Lno=Length[NonRelOrbitals];
	NonRelConfigs={};
	RelConfigs={};
	For[k=1,k<=Length[Shells],k=k+1,
	{
		list=ShellList[Shells[[k]]];
		NonRelOrbits=ConstantArray[0,Lno];
		RelOrbits=ConstantArray[0,Lro];
		For[i=1,i<=2,i=i+1,
		{
			n=Count[list,RelOrbitals[[i]]];
			NonRelOrbits[[i]]=n;
			RelOrbits[[i]]=n;
		}];
		j=3;
		For[i=3,i<=Lro,i=i+2,
		{
			n1=Count[list,RelOrbitals[[i]]];
			n2=Count[list,RelOrbitals[[i+1]]];
			NonRelOrbits[[j]]=n1+n2;
			RelOrbits[[i]]=n1;
			RelOrbits[[i+1]]=n2;
			j=j+1;
		}];
		AppendTo[NonRelConfigs,NonRelOrbits];
		AppendTo[RelConfigs,RelOrbits];
	}];
	CONFIGS=DeleteDuplicates[NonRelConfigs];
	NumPROJECT=ConstantArray[0,Length[CONFIGS]];
	NumCONFIGS=ConstantArray[0,Length[CONFIGS]];
	For[i=1,i<=Length[CONFIGS],i=i+1,
	{
		pos=Position[NonRelConfigs,CONFIGS[[i]]];
		NumPROJECT[[i]]=Total[Num[[pos\[Transpose][[1]]]]];
		NumCONFIGS[[i]]=Length[pos];
	}];
};
Return[{CONFIGS,NumPROJECT,NonRelConfigs,NumCONFIGS,RelConfigs}]
]


ConigurationAverageJ[CONFIGS_]:=Module[{AverageJ,L,m,config,shellj,k,g,ml,J,jcount,i,j,a,Aa,Ba,Ca,sum},
{
	AverageJ=ConstantArray[0,Length[CONFIGS]];
	L=Length[NonRelOrbitals];
	For [m=1,m<= Length[CONFIGS],m=m+1,
	{
		config=CONFIGS[[m]];
		shellj={};
		For[k=1,k<=L,k=k+1,
		{
			g=config[[k]];
			If[g==0||g==MaxElectrons[[k]],Continue[]];
			ml=Range[-l[[k]],l[[k]]];
			J=ConstantArray[0,Length[ml]*Length[spin]];
			jcount=1;
			For[i=1,i<=Length[ml],i++,
			{
				For[j=1,j<=Length[spin],j++,
				{
					a={ml[[i]],spin[[j]]};
					J[[jcount]]=a;
					jcount++1;
				}]
			}];
			Aa=Subsets[J,{g}];
			Ba=Total[Aa,{3}];
			Ca=Total[Ba,{2}];
			AppendTo[shellj,Ca];
		}];
		sum=Abs[Total[Tuples[shellj],{2}]];
		AverageJ[[m]]=Mean[sum];
	}];
};
Return[AverageJ]
];


(*Absorption Spectrum \[Rule] Opposite of emission spectrum*)
(* initial \[TwoWayRule] final *)
DiscreteLevels[FILE_,Temp_]:=Block[{Emin,Four,Two,CIinitial,CIfinal,CIindi,CIindf,j,EMAX,J,En,Einstein,CI,S,L,Egap,i,T,Jinitial,Jfinal,Pinitial,Pfinal,Indinitial,Indfinal,Jposi,Jposf,Einitial,EINITIAL,Efinal,EFINAL,Intensity,Z,pos,Boltz,Jupp,Jlow,Eupp,Elow,g,A},
{
	CI=ReadCI[FILE];
	{Jinitial,Pinitial,CIinitial,Jfinal,Pfinal,CIfinal,S}=ReadTransitions[FILE];

	L=Length[Jinitial];
	
	En=CI[[All,All,4,All,2]];
	En=En eV;
	Emin=Min[En];
	En=En-Emin;
	J=CI[[All,All,1]];
	Z=Total[Flatten[(2J+1)Exp[-En/Temp]]];
	

	Jposi=Floor[Jinitial]+1;
	Jposf=Floor[Jfinal]+1;
	CIindi=CIinitial+1;
	CIindf=CIfinal+1;
	Four=ConstantArray[4,Length[Jinitial]];
	Two=ConstantArray[2,Length[Jinitial]];
	EINITIAL=Table[CI[[Pinitial[[i]],Jposi[[i]],Four[[i]],CIindi[[i]],Two[[i]]]],{i,1,L}];
	EFINAL=Table[CI[[Pfinal[[i]],Jposf[[i]],Four[[i]],CIindf[[i]],Two[[i]]]],{i,1,L}];
	Egap=Abs[EINITIAL-EFINAL]eV;
	EMAX=Map[Max,Flatten[{EINITIAL,EFINAL},{{2},{1}}],{1}];
	
	A=(2.0261*10^18)(*/g*)*(Egap)^3/(4.1357*10^-15*3*10^8*10^10)^3 S;
	Einstein=A(**g*);
	Boltz=(*g*) Exp[-(EMAX eV-Emin)/Temp];
	Intensity=(Egap)*A/(4\[Pi] Z)*Boltz;

};
Return[{Egap,S,EINITIAL,EFINAL,Intensity,Einstein}]
]

DiscreteSpectrum[Energy_,S_,w_,G_]:=Module[{Normalise,L,Spec,Spectrum,i,j},
{
	L=Length[S];
	Normalise=1/(\[Pi] G);
	(*Print["Convolution Started"];*)
	Spectrum=Sum[Normalise*G^2/((w-Energy[[i]])^2+G^2)*S[[i]],{i,1,L},Method->"Procedural"];
	(*Print["Convolution Finished"];*)
};
Return[Spectrum]
];



ComposeSingle[STR_,Orbitals_,EnergySingle_]:=Module[{OmegaFinal,OmegaFinalPos,OmegaInitialPos,OmegaInitial,Dmel,DeltaOmega,ShellInitial,ShellFinal,jInitial,nInitialPos,jFinal,nFinalPos},
{
	Dmel=ToExpression[STR[[5]]]^2//N;
	ShellInitial=STR[[1]];
	ShellFinal=STR[[3]];
	OmegaInitialPos=Position[Orbitals,ShellInitial][[1]];
	OmegaInitial=EnergySingle[[OmegaInitialPos[[1]]]];
	OmegaFinalPos=Position[Orbitals,ShellFinal][[1]];
	OmegaFinal=EnergySingle[[OmegaFinalPos[[1]]]];
	DeltaOmega=OmegaInitial-OmegaFinal;

	jInitial=Jsingle[{ShellInitial}];
	nInitialPos=Position[RelOrbitals,ShellInitial][[1]];

	jFinal=Jsingle[{ShellFinal}];
	nFinalPos=Position[RelOrbitals,ShellFinal][[1]];
};
Return[{Dmel,DeltaOmega,ShellInitial,ShellFinal,jInitial,nInitialPos,jFinal,nFinalPos}];
];

SingleTransitionData[FILE_]:=Module[{Orbitals,EnergySingle,Transitions,T,Col1,Col2,Col3,Col4,Col5,TransitionsInverse,TransitionsBoth,Ltran},
{
	{Orbitals,EnergySingle}=ReadOrbitals[FILE];
	EnergySingle=EnergySingle eV;
	Transitions=ReadSingleTransition[FILE];
	Col1=Transitions\[Transpose][[1]];
	Col2=Transitions\[Transpose][[2]];
	Col3=Transitions\[Transpose][[3]];
	Col4=Transitions\[Transpose][[4]];
	Col5=Transitions\[Transpose][[5]];
	TransitionsInverse={Col3,Col2,Col1,Col4,Col5}\[Transpose];
	TransitionsBoth=Join[Transitions,TransitionsInverse];
	Ltran=Length[TransitionsBoth];
	T=Table[ComposeSingle[TransitionsBoth [[i]],Orbitals,EnergySingle],{i,1,Ltran}];
};
Return[T];
]


E1Transitions[i_,j_,Shells_,T_,nAverage_]:=Module[{E1,A,B,Diff1,Diff2,exp,p1,p2,Int,n1,n2},
{
	If[i==j,
	{
		E1=0;
		exp=0;
	},{
		A=ShellList[Shells[[i]]];
		B=ShellList[Shells[[j]]];
		{Diff1,Diff2}=ShellDifferences[A,B];
		If[Length[Diff1]>1||Length[Diff2]>1,
		{
			E1=0;
			exp=0
		},{
			p1=Position[T\[Transpose][[3]],Diff1[[1]]];
			p2=Position[T\[Transpose][[4]],Diff2[[1]]];
			Int=Intersection[p1,p2];
			If[Int=={},
			{
				E1=0;
				exp=0;
			},{
				E1=Int[[1]][[1]];
				n1=nAverage[[i]][[T[[E1]][[6]]]][[1]];
				n2=nAverage[[i]][[T[[E1]][[8]]]][[1]];
				exp=n1(1-n2)//N;
			}];
		}]
	}];
};
Return[{E1,exp}]
];

ExpVal[k_,Ck2_,E1_]:=Module[{EXP,P,P2,exp,CIstate,w},
{
	P=Position[E1,k];
	If[P=={},
	{
		EXP=0;
	},{
		P2=P;
		P2[[All,3]]=2;
		exp=Extract[E1,P2];
		CIstate=P[[All,1]];
		w=Ck2[[CIstate]];
		EXP=Total[w*exp];
	}];
};
Return[EXP]
]

ListIndicies[L_]:=Module[{TotalLines,iPsi,jPsi,ind,i,j},
{
	TotalLines=L (L-1)/2;
	iPsi=ConstantArray[0,TotalLines];
	jPsi=ConstantArray[0,TotalLines];
	ind=0;
	For[i=1,i<=L,i++,
	{
		For[j=i+1,j<=L,j++,
		{
			ind++;
			iPsi[[ind]]=i;
			jPsi[[ind]]=j;
		}];
	}];
};
Return[{iPsi,jPsi}];
]


Main[ind_,Cutoff_,Renorm_,G_,T_,Light_,indicies_,Dspacing2_,Jin_,NumLines_,E1_,Ltran_,Components2_]:=Module[{jPsi,iPsi,Ck2,Strength,Lorenzian,EXP,SUM},
{
	iPsi=indicies[[1]][[ind]];
	jPsi=indicies[[2]][[ind]];

	(*Lorenzian=1/(2\[Pi])*(G+G)/((T\[Transpose][[2]]-Light[[ind]])^2+(G+G)^2/4);*)
	Lorenzian=1/(2\[Pi])*(G+G)/((T\[Transpose][[2]]-Light[[ind]])^2+(G+G)^2/4)/Renorm*HeavisideTheta[Cutoff-(T\[Transpose][[2]]-Light[[ind]])]*HeavisideTheta[Cutoff+(T\[Transpose][[2]]-Light[[ind]])];
	(*Gaussian=1/Sqrt[Sigma^2*2\[Pi]]*Exp[- (T\[Transpose][[2]]-Light[[ind]])^2/(2Sigma^2)];*)

	Ck2=Components2[[jPsi]];
	EXP=Table[ExpVal[k,Ck2,E1],{k,1,Ltran}];
	SUM=Total[T\[Transpose][[1]]*Lorenzian*EXP];
	(*SUM=Total[T\[Transpose][[1]]*Gaussian*EXP];*)

	Strength=SUM*Dspacing2[[ind]]*(2*Jin[[ind]]+1)*NumLines[[ind]]/3;

};
Return[Strength]
]



End[];
EndPackage[];
