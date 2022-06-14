(* ::Package:: *)

(*
This program is built based on the technique presented in the following articles:

[1] Z. Zhang, "A Flexible New Technique for Camera Calibration", IEEE Transactions on Pattern Analysis and Machine Intelligence, 2000.
[2] W. Burger, "Zhang's Camera Calibration Algorithm: In-Depth Tutorial and Implementation", University of Applied Sciences Upper Austria, HGB16-05, 2016.

*)


funcZhang[pAs_,pBs_,LL_]:=Module[{norm,T,getZeroSpace,sg,twolines,MM,hh,hh1a,hh2a,hh3a,hh1b,hh2b,hh3b,h2v,VV,B11,B12,B13,B22,B23,B33,v0,lmd,alp,bet,gmm,u0,AA,TTA,TTB,TTAB,funcHtoT},
(*
	This function accepts the information on flat patterns A and B and provides the position/attitude relation between them. 
	The patterns A and B are 4x4 grids and define the coordinate frames A and B, respectively. 
	The arguments pAs and pBs are 4x4x2 arrays, representing the patterns A and B, respectively.
	The argument LL is the physical length of the grid pitch of the patterns A and B.
	The outputs are:
	  TTAB: 4x4 matrix: Frame B seen from Frame A
	  AA: 3x4 matrix: Camera Frame seen from Image
	  TTA: 4x4 matrix: Frame A seen from Camera Frame
	  TTB: 4x4 matrix: Frame B seen from Camera Frame
*)
	sg=-1; (* To make the camera frame to be right-handed. *)
	norm[a_]:=Sqrt[a.a];
	T=Transpose;
	getZeroSpace[MM_]:=Module[{absMax,ss,hhh},
		absMax[a_]:=Max[Table[Abs[a[[i]]],{i,Length[a]}]];
		ss=DiagonalMatrix[Table[1/absMax[T[MM][[i]]],{i,Length[T[MM]]}]];
		hhh=T[SingularValueDecomposition[MM.ss][[3]]][[-1]];
		(*SingularValueList[MM.ss];//Print*)
		Return[ss.hhh];
	];
	h2v[hi_,hj_]:={hi[[1]]*hj[[1]],hi[[1]]*hj[[2]]+hi[[2]]*hj[[1]],hi[[2]]*hj[[2]],hi[[3]]*hj[[1]]+hi[[1]]*hj[[3]],hi[[3]]*hj[[2]]+hi[[2]]*hj[[3]],hi[[3]]*hj[[3]]};
	twolines[XYk_,uvk_]:={
		{-XYk[[1]],-XYk[[2]],-1,0,0,0,sg*uvk[[1]]*XYk[[1]],sg*uvk[[1]]*XYk[[2]],sg*uvk[[1]]},
		{0,0,0,-XYk[[1]],-XYk[[2]],-1,sg*uvk[[2]]*XYk[[1]],sg*uvk[[2]]*XYk[[2]],sg*uvk[[2]]}};
	MM=ArrayFlatten[Flatten[Table[Table[twolines[{{-LL,0,LL,2LL}[[-j]],{-LL,0,LL,2LL}[[-i]]},pAs[[i,j]]],{i,1,4}],{j,1,4}],2]];
	hh=getZeroSpace[MM];
	hh1a=hh[[{1,4,7}]];
	hh2a=hh[[{2,5,8}]];
	hh3a=hh[[{3,6,9}]];
	MM=ArrayFlatten[Flatten[Table[Table[twolines[{{-LL,0,LL,2LL}[[-j]],{-LL,0,LL,2LL}[[-i]]},pBs[[i,j]]],{i,1,4}],{j,1,4}],2]];
	hh=getZeroSpace[MM];
	hh1b=hh[[{1,4,7}]];
	hh2b=hh[[{2,5,8}]];
	hh3b=hh[[{3,6,9}]];
	VV={
		h2v[hh1a,hh2a],h2v[hh1a,hh1a]-h2v[hh2a,hh2a],
		h2v[hh1b,hh2b],h2v[hh1b,hh1b]-h2v[hh2b,hh2b],
		{0,1,0,0,0,0},{1,0,-1,0,0,0}};
	{B11,B12,B22,B13,B23,B33}=getZeroSpace[VV];
	B13=sg* B13;(*Si son negativos los 3 salen imaginarios*)
	B23=sg* B23;
	B33=Abs[B33];
	v0=N[(B12*B13-B11*B23)/(B11*B22-B12^2)];
	lmd=N[B33-(B13^2+v0*(B12*B13-B11*B23))/B11];
	alp=N[Sqrt[lmd/B11]];
	bet=N[Sqrt[lmd*B11/(B11*B22-B12^2)]];
	gmm=N[-((B12*alp^2*bet)/lmd)];
	u0=N[(gmm*v0)/alp-(B13*alp^2)/lmd];
	AA={{alp,gmm,u0},{0,bet,v0},{0,0,sg}};
	funcHtoT[AA_,hh1_,hh2_,hh3_]:=Module[{iAA,rAx,ss,rAy,rAz,rAt},
		iAA=Inverse[AA];
		rAx=iAA.hh1;
		ss=1./norm[rAx];
		rAx=ss*rAx;
		rAy=ss*iAA.hh2;
		rAz=Cross[rAx,rAy];
		rAt=ss*iAA.hh3;
		Return[T[{rAx,rAy,rAz,rAt}]];
	];
	TTA=funcHtoT[AA,hh1a,hh2a,hh3a];
	TTB=funcHtoT[AA,hh1b,hh2b,hh3b];
	AA=ArrayFlatten[{{AA,{{0},{0},{0}}}}];
	TTA=ArrayFlatten[{{TTA},{{{0,0,0,1}}}}];
	TTB=ArrayFlatten[{{TTB},{{{0,0,0,1}}}}];
	TTAB=Inverse[TTA].TTB;
	Return[{TTAB,AA,TTA,TTB}];
];
