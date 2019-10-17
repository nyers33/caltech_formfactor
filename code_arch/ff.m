(*
	Copyright (c) 1992, 1993 Peter Schr\"oder
		Department of Computer Science
		35 Olden St.
		Princeton University
		Princeton, NJ 08544-2087
		ps@cs.princeton.edu

	Peter Schr\"oder hereby grants anyone permission to use, copy and
modify this program, so long as this notice in its entirety is attached.
This permission does not include the right to redistribute or sublicense
for consideration, or commercially exploit in any way, any modification of
this program, any derivative work made from this program, or any
incorporation of this program in another program, without the express
written consent of Peter Schr\"oder.
                          
	PETER SCHR\"ODER PROVIDES ABSOLUTELY NO WARRANTY OF ANY KIND WITH
RESPECT TO THIS PROGRAM, INCLUDING WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE.  THE ENTIRE RISK AS TO THE QUALITY AND
PERFORMANCE OF THIS PROGRAM IS WITH THE USER.  IN NO EVENT WILL PETER
SCHR\"ODER BE LIABLE TO ANYONE FOR DAMAGES ARISING OUT OF THE USE OF THIS
PROGRAM, INCLUDING, WITHOUT LIMITATION, DAMAGES RESULTING FROM LOST DATA OR
LOST PROFITS, OR FOR ANY SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES.
*)

Remove[pair,bilcoeff,fourc,integralplanar,G,H,
	l1,l2,l3,l4,logs,lcis,logselect,kpart,dilog,logpart,
	logint,integral,mag,area,cross,FormFactor]

pair[l_List]:=
	Block[{pi,pj,ei,ej,di,dj},

		(*
		   given two lines compute c_{0,1,2,3,4,5}
		   Mathematica arrays are offset 1 so they
		   are actually called c[[1,2,3,4,5,6]]...
		*)
		pj=l[[1,1]];
		pi=l[[2,1]];
		ej=l[[1,2]]-l[[1,1]];
		ei=l[[2,2]]-l[[2,1]];
		dj=ej/mag[ej];
		di=ei/mag[ei];
		Simplify[{mag[ej],-2 di.dj,mag[ei],-2 dj.(pi-pj),
			2 di.(pi-pj),(pi-pj).(pi-pj)}]
	]

bilcoeff[c_List]:=
	Block[{const,theta,phi,argdiff,eps=10^-6},

		(*
		   determine whether the biquadratic form
		   factors into the product of two bilinear
		   forms. This corresponds to the case that
		   the two edges lie in the same plane
		   Since we always assume the coefficients
		   of s^2 and t^2 to be 1 we may write the
		   product of the bilinear forms as

		   (s Exp[I theta]+t Exp[I phi]+const)
		   (s Exp[-I theta]+t Exp[-I phi]+const)

		   From this follows that Sqrt[c[[6]]]
		   must be const
		*)
		const=Sqrt[c[[6]]];
		(*
		   Furthermore, if c==0 then c[[4]] and
		   c[[5]] must also be zero! In this
		   case we are underconstrained. What's more
		   Abs[c[[2]]]<=2 since

		   ArcCos[c[[2]]/2]=theta-phi

		   since both c[[4]] and c[[5]] are zero
		   we arbitrarily pick theta==0
		*)
		argdiff=ArcCos[c[[2]]/2];
		If[N[Abs[const]]<eps,
			theta=0;
			phi=-argdiff,
			(*
			   now we have to satisfy
			*)
			theta=ArcCos[c[[4]]/(2const)];
			phi=ArcCos[c[[5]]/(2const)]
		];
		(*
		   we still have the dofs to contend with
		   which arise from the fact that the
		   ArcCos is only determined up to sign
		*)
		If[N[Abs[Arg[Exp[I (theta-phi)]]-argdiff]]<eps,
			{theta,phi,const},
		If[N[Abs[Arg[Exp[I (theta+phi)]]-argdiff]]<eps,
			{theta,-phi,const},
		If[N[Abs[Arg[Exp[I (-theta-phi)]]-argdiff]]<eps,
			{-theta,phi,const},
		If[N[Abs[Arg[Exp[I (-theta+phi)]]-argdiff]]<eps,
			{-theta,-phi,const},
			{}
			]
			]
			]
		]
	]

(*
   evaluate f at the four corners
*)
fourc[f_,s_,t_]:=f[s,t]-f[s,0]-f[0,t]+f[0,0]

integralplanar[bil_List,c0_,c2_]:=
	Block[{h,int},

		(*
		   if the edges share a plane this will
		   be the entire answer
		*)
		h[s_,t_]=s Exp[I bil[[1]]]+t Exp[I bil[[2]]]+bil[[3]];
		(*
		   if h[s,t] goes to zero at any one of the corners
		   of integration, that corner will be zero (a limit
		   argument that, unfortunately, Mathematica is not
		   capable of)
		*)
		int[s_,t_]=If[N[h[s,t]==0],
			0,
			Re[h[s,t]^2/Exp[I (bil[[1]]+bil[[2]])]*
				(Log[h[s,t]]-3/2)]
		];
		int[c0,c2]-int[c0,0]-int[0,c2]+int[0,0]
	]

G[a_,b_,c_,t_]:=
	Block[{q,qp,d=Sqrt[4a c - b^2]},

		(*
		   the integral of `log[a t^2+b t+c]'
		*)
		q[y_]=a y^2 + b y + c;
		qp[y_]=2a y+b;
		qp[t]/(2a) Log[q[t]]-2t+d/a*ArcTan[qp[t]/d]
	]

H[a_,b_,c_,t_]:=
	Block[{q,qp,d=Sqrt[4a c - b^2]},

		(*
		   the integral of `t log[a t^2+b t+c]'
		*)
		q[y_]=a y^2 + b y + c;
		qp[y_]=2a y+b;
		(2a^2 t^2+2a c-b^2)/(4a^2) Log[q[t]]-
			t(a t-b)/(2a) - b d/(2a^2) ArcTan[qp[t]/d]
	]

(*
   the four different continuous extensions of the logarithm
   that we'll need depending on the integration path
*)
l1[x_]:=Log[x]

l2[x_]:=If[Im[x]<=0&&Re[x]<=0,Log[x]+2I Pi,Log[x]]

l3[x_]:=If[Im[x]>=0,Log[x]-2I Pi,Log[x]]

l4[x_]:=If[Im[x]>=0&&Re[x]<=0,Log[x]-2I Pi,Log[x]]

logs={l1,l2,l3,l4}

lcis[a_,b_]:=
	Block[{xdisc=Re[b]^2-Abs[b]^2+Abs[a]^2,
		ydisc=Im[b]^2-Abs[b]^2+Abs[a]^2,
		x1,x2,y1,y2},

		(*
		   intersect the real and the imaginary axis
		   with a circle of radius a, centered at b
		*)
		x1=Re[b]+Sqrt[xdisc];
		x2=Re[b]-Sqrt[xdisc];
		y1=Im[b]+Sqrt[ydisc];
		y2=Im[b]-Sqrt[ydisc];

		(*
		   return the solutions, if any, in the order
		   real axis, imaginary axis
		*)
		{If[xdisc>=0,{x1,x2},{}],If[ydisc>=0,{y1,y2},{}]}
	]

logselect[a_,b_,l_,u_]:=
	Block[{maxpath,minpath,posreal,negreal,posimg,negimg,
		seghit,sols},
		(*
		   a is the radius
		   b is the center
		   l is lower
		   u is upper
		   any argument has the form
		      z=b+a*Exp[I(Arg[l]+t(Arg[u]-Arg[l]))]
		   resolve the multivaluedness accordingly

		   This is your basic graphics problem. We
		   have a unit circle centered at b and scaled by a and ask
		   whether the real axis or the imaginary axis intersects it,
		   and if so, in the segment of the circle that
		   we will traverse during the integration.
		*)
		maxpath=Max[Arg[l],Arg[u]];
		minpath=Min[Arg[l],Arg[u]];

		seghit[z_]:=minpath<=Arg[(z-b)/a]<=maxpath;

		sols=lcis[a,b];

		(* check for negative real axis *)
		negreal=(Length[sols[[1]]]==2)&&
			((Re[sols[[1,1]]]<=0 && seghit[sols[[1,1]]])||
			(Re[sols[[1,2]]]<=0 && seghit[sols[[1,2]]]));
		(* check for positive real axis *)
		posreal=(Length[sols[[1]]]==2)&&
			((Re[sols[[1,1]]]>=0 && seghit[sols[[1,1]]])||
			(Re[sols[[1,2]]]>=0 && seghit[sols[[1,2]]]));
		(* check for negative imaginary axis *)
		negimg=(Length[sols[[2]]]==2)&&
			((Re[sols[[2,1]]]<=0 && seghit[I sols[[2,1]]])||
			(Re[sols[[2,2]]]<=0 && seghit[I sols[[2,2]]]));
		(* check for positive imaginary axis *)
		posimg=(Length[sols[[2]]]==2)&&
			((Re[sols[[2,1]]]>=0 && seghit[I sols[[2,1]]])||
			(Re[sols[[2,2]]]>=0 && seghit[I sols[[2,2]]]));

		If[!negreal,
			1,
			If[!negimg,
				2,
				If[!posimg,
					4,
					If[!posreal,
						3,
					]
				]
			]
		]
	]

kpart[k_,l_,u_]:=
	Block[{int},

		(* int[t]=Integrate[t^2/(1-t^2)^3,t] *)
		int[t_]:=t/(4(t^2-1)^2)+t/(8(t^2-1))+Log[(t-1)/(t+1)]/16;

		I(2k+1)Pi(int[u]-int[l])
	]

dilog[c_,r_,l_,u_]:=
	Block[{sols,maxpath,minpath,seghit,cross,p,parg,towardsu,towardsl,
		pu,pl,start,end},

		(*
		   return the points where the argument of the
		   log is on the real axis
		*)
		sols=lcis[r,c][[1]];

		(*
		   since we don't know the ordering
		*)
		maxpath=Max[Arg[l],Arg[u]];
		minpath=Min[Arg[l],Arg[u]];

		(*
		   given such a point on the real axis
		   where does it lie on the unit cirle
		   (of which we are interested in the
		   segment between l and u)
		*)
		seghit[z_]:=minpath<=Arg[(z-c)/r]<=maxpath;

		(*
		   only if there are any solutions,
		   the solutions are on the negative real axis,
		   they are within the segment of interest,
		   there is only one intersection (we actually
		   cross)
		*)
		cross=(Length[sols]==2)&&
			Xor[Re[sols[[1]]]<=0 && seghit[sols[[1]]],
			Re[sols[[2]]]<=0 && seghit[sols[[2]]]];

		(*
		   for the dilog we have an argument of the form
		   1-z, compute it
		*)
		p=If[Re[sols[[1]]]<=0 && seghit[sols[[1]]],
			1-sols[[1]],
			1-sols[[2]]
		];

		(*
		   put this argument on the unit circle
		   remembering that we already know that it
		   is within the contour over which we are
		   integrating
		*)
		parg=Arg[((1-p)-c)/r];

		(*
		   we want to cut the path at that point
		*)
		towardsu=If[Arg[u]<parg,parg-10^-10,parg+10^-10];
		towardsl=If[Arg[l]<parg,parg-10^-10,parg+10^-10];

		(*
		   now that we have these cut points, put them
		   back into the coordinate frame we'll use
		   for the dilog, remembering that the unit circle
		   segment between l and u will be subjected to
		   the transformation 1-(c+r t)
		*)
		pu=1-(c+r Exp[I towardsu]);
		pl=1-(c+r Exp[I towardsl]);
		start=1-(c+r l);
		end=1-(c+r u);

		(*
		   which particular extension of the standard
		   log are we integrating over?
		*)
		k=logselect[r,c,l,u];

		PolyLog[2,end]-PolyLog[2,start]+
		If[cross,
			(*
			   correction factor based on the particular
			   extension of the log we chose
			*)
			If[k==2,
				If[Im[pl]>0,
					-2Pi I(2 Log[p]-Log[start]),
					-2Pi I Log[end]
				],
				If[k==4,
					If[Im[pu]<0,
						-2Pi I(2 Log[p]-Log[end]),
						-2Pi I Log[start]
					],
					0
				]
			],
			0
		]
	]

logpart[b_,l_,u_]:=
	Block[{p1,p2,p3,int},

		(*
		   this integrates `t^2/(1-t^2)^3 Log[t+b]'
		   between l and u
		*)
		p1=logs[[logselect[1,b,l,u]]];
		p2=logs[[logselect[-1/(1+b),1/(1+b),l,u]]];
		p3=logs[[logselect[1/(1-b),1/(1-b),l,u]]];
		int[t_]:=-p1[t-1]b/(1+b)^2-p1[1+t]b/(b-1)^2+
			p1[b+t]2(b+t)(1+b t)((b-t)^2+(b t-1)^2)/
			((b^2-1)^2(t^2-1)^2)+2(b-t)/((b^2-1)(t^2-1))+
			p2[(1-t)/(1+b)]p1[b+t]-p3[(1+t)/(1-b)]p1[b+t];

		1/16(int[u]-int[l]+
			dilog[1/(1+b),-1/(1+b),l,u]-
			dilog[1/(1-b),1/(1-b),l,u])
	]

logint[s_,c15_,c16_,c17_,c18_,l_,u_]:=
	Block[{k,lan,p1,p2,p3,p4,tan,check,eps=10^-6},
			
		(*
		   we need to determine the value of k
		   when going from the arctan to the log
		   representation. We simply guess and stick
		   in two points to see what we get
		*)
		p1=logs[[logselect[1,-c17[s],l,u]]];
		p2=logs[[logselect[1,-c18[s],l,u]]];
		p3=logs[[logselect[1,c17[s],l,u]]];
		p4=logs[[logselect[1,c18[s],l,u]]];
		(*
		   the integrand when written with logs
		*)
		lan[t_,k_]:=t^2/(1-t^2)^3/(2I)*(I(2k+1)Pi+
			p1[t-c17[s]]+p2[t-c18[s]]-
			p3[t+c17[s]]-p4[t+c18[s]]);
		(*
		   the original integrand
		*)
		tan[t_]=t^2/(1-t^2)^3*
			ArcTan[(Conjugate[c16[s]] t^2-c16[s])/(c15 t)];
		check[k_]:=N[Abs[tan[u]-lan[u,k]]+Abs[tan[l]-lan[l,k]]];
		k=If[check[0]<eps,
			0,
			If[check[1]<eps,
				1,
				If[check[-1]<eps,
					-1,
					(* oh, oh... *)
					0
				]
			]
		];

		kpart[k,l,u]+
			logpart[-c17[s],l,u]+
			logpart[-c18[s],l,u]-
			logpart[c17[s],l,u]-
			logpart[c18[s],l,u]
	]

integral[c_List]:=
	Block[{bc,c6,c7,r,c10,c11,c12,c13,c14,c15,c16,d,c17,c18,
		firstpart,upper,lower},

		(*
		   see whether it can be factored into two bilinear
		   forms
		*)
		bc=bilcoeff[c];

		If[Length[bc] != 0,
			integralplanar[bc,c[[1]],c[[3]]],

			(*
			   the general case
			   compute all the coefficients from the paper
			*)
			c6[t_]=c[[2]] t + c[[4]];
			c7[t_]=t^2 + c[[5]] t + c[[6]];
			r[t_]=4 c7[t]-c6[t]^2;
			c10=Coefficient[r[t],t,2];
			c11=Coefficient[r[t],t,1];
			c12=Coefficient[r[t],t,0];
			c13=(c11 - Sqrt[c11^2-4c10 c12])/(2c10);
			c14=Sqrt[c11^2-4c10 c12]/c10;
			c15=Sqrt[c10] c14;
			c16[s_]=c[[2]] c13-c[[4]]-2 s;
			(*
			   new limits of integration
			   for sanity we'll make them both situated
			   in the upper half complex plane.
			*)
			upper=Sqrt[(c13+c[[3]])/(Conjugate[c13]+c[[3]])];
			upper=If[Im[N[upper]]<0,-upper,upper];
			lower=Sqrt[c13/Conjugate[c13]];
			lower=If[Im[N[lower]]<0,-lower,lower];
			c17[s_]=(-c15+
				Sqrt[c15^2-4 c16[s] Conjugate[c16[s]]])/
				(-2Conjugate[I c16[s]]);
			c18[s_]=(-c15-
				Sqrt[c15^2-4 c16[s] Conjugate[c16[s]]])/
				(-2Conjugate[I c16[s]]);
			(* equation 5 *)
			firstpart[s_,t_]=
				(s+c[[4]]/2)*
				G[1,c[[5]]+c[[2]] s,c[[6]]+c[[4]] s + s^2,t]+
				c[[2]]/2*
				H[1,c[[5]]+c[[2]] s,c[[6]]+c[[4]] s + s^2,t];

			fourc[firstpart,c[[1]],c[[3]]]-2c[[1]]c[[3]]-
			I c14 c15(
			(*
			   often c16[] is zero at one of the extremes
			   indicating that the
			   whole inverse tangent integral is zero.
			   Mathematica doesn't know this and dutifully
			   finds divisions by zero, etc. Hence we test
			   for it explicitely
			*)
			If[c16[c[[1]]]!=0,
				logint[c[[1]],c15,c16,c17,c18,lower,upper],
				0]-
			If[c16[0]!=0,
				logint[0,c15,c16,c17,c18,lower,upper],
				0
			])
		]
	]

mag[v1_]:=Sqrt[v1.v1]

cross[a_,b_]:=
	{a[[2]]*b[[3]]-a[[3]]*b[[2]],
	 a[[3]]*b[[1]]-a[[1]]*b[[3]],
	 a[[1]]*b[[2]]-a[[2]]*b[[1]]}

area[p_]:=
	Block[{i,ni,d1={0,0,0},d2={0,0,0},a={0,0,0}},

		ni[k_]:=Mod[k,Length[p]]+1;
		For[i=1,i<Length[p],i+=1,
			d1=p[[ni[i]]]-p[[1]];
			d2=p[[ni[ni[i]]]]-p[[1]];
			a+=cross[d1,d2]
		];
		
		mag[a]/2
	]

FormFactor[p1_List,p2_List]:=
	Block[{ni,nj,i,j,ff=0,c={},eps=10^-6},

		ni[k_]:=Mod[k,Length[p1]]+1;
		nj[k_]:=Mod[k,Length[p2]]+1;
		For[i=1,i<=Length[p1],i+=1,
			For[j=1,j<=Length[p2],j+=1,
				c=pair[{{p1[[i]],p1[[ni[i]]]},
					{p2[[j]],p2[[nj[j]]]}}];
				If[Abs[c[[2]]]>eps,
					ff-=c[[2]] integral[c],
				]
			]
		];

		ff/(8 Pi area[p1])
	]
