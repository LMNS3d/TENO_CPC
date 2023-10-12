//-------------------------------------------------------------------------------------------------
//		Subroutine for calculating the cell-interface flux in x direction
//-------------------------------------------------------------------------------------------------
void CPUFluid::CharacteristicF_x(int i, int j, int k, Real UI[Emax][Xmax][Ymax][Zmax])
{
	//Scalar u, f, f_plus, f_minus
	Real uf[10], ff[10], pp[10], mm[10], f_flux, _p[Emax][Emax];
	fh1=0.0; fh2=0.0; fh3=0.0; fh4=0.0; fh5=0.0;
	rhom=0.0; uui=0.0; vvi=0.0;wwi=0.0;ppi=0.0;enti=0.0;
	uuip=0.0; vvip=0.0; wwip=0.0; ppip=0.0;entip=0.0;
	ft1=0.0;ft2=0.0;ft3=0.0;ft4=0.0; ft5=0.0; ft6=0.0;
	uvs1=0.0;uvs2=0.0;uvs3=0.0; uvs4=0.0; uvs5=0.0; uvs6=0.0; uv_part=0.0;
	Real xishu[4]={0.8441336457279522, -0.2441336457279522, 0.0570096576929319, -0.006723831837710873};
	// start
	int indicator_x = 0;
	for (int m = 0; m < stencil_size; m++) 
		uf[m] = UI[0][m + i - stencil_P][j][k];

	indicator_x = AMAX1(indicator_x, indicator(&uf[stencil_P]));

	if (DIM_X + DIM_Y + DIM_Z == 1)
	{
		for (int m = 0; m < stencil_size; m++)
			uf[m] = UI[Emax - 1][m + i - stencil_P][j][k];

		indicator_x = AMAX1(indicator_x, indicator(&uf[stencil_P]));
	}
	// for statistics

	fcells[i][j][k]->tag_linear_x = indicator_x;

	// end the statistics

	if (indicator_x == 0)
	{
		// return the skew-symmetric-splitting flux
		F_x_wall[0][i][j][k]=0.0;
		F_x_wall[1][i][j][k]=0.0;
		F_x_wall[2][i][j][k]=0.0;
		F_x_wall[3][i][j][k]=0.0;
		F_x_wall[4][i][j][k]=0.0;
		ft1 = 0.0;
		ft2 = 0.0;
		ft3 = 0.0;
		ft4 = 0.0;
		ft5 = 0.0;
		ft6 = 0.0;
		for(int l=1; l<5; l++){
			uvs1 = 0.0;
			uvs2 = 0.0;
			uvs3 = 0.0;
			uvs4 = 0.0;
			uvs5 = 0.0;
			uvs6 = 0.0;
			for(int m=0; m<l; m++){
				rhom  = rho[i-m][j][k] + rho[i-m+l][j][k];

				uui   = u[i-m][j][k];
				vvi   = v[i-m][j][k];
				wwi   = w[i-m][j][k];
				ppi   = p[i-m][j][k];
				enti  = U[4][i-m][j][k];

				uuip  = u[i-m+l][j][k];
				vvip  = v[i-m+l][j][k];
				wwip  = w[i-m+l][j][k];
				ppip  = p[i-m+l][j][k];
				entip = U[4][i-m+l][j][k];

				uv_part = (uui+uuip) * rhom;
				uvs1 = uvs1 + uv_part * (2.0);
				uvs2 = uvs2 + uv_part * (uui+uuip)+(ppi+ppip)* (4.0);
				uvs3 = uvs3 + uv_part * (vvi+vvip);
				uvs4 = uvs4 + uv_part * (wwi+wwip);
				uvs5 = uvs5 + (uui+uuip) * (enti+entip)* (2.0)+ (uui+uuip) * (ppi+ppip)* (2.0);
				
			}
			ft1 = ft1 + xishu[l-1]*uvs1;
			ft2 = ft2 + xishu[l-1]*uvs2;
			ft3 = ft3 + xishu[l-1]*uvs3;
			ft4 = ft4 + xishu[l-1]*uvs4;
			ft5 = ft5 + xishu[l-1]*uvs5;
		}
		fh1 = 0.25*ft1;
		fh2 = 0.25*ft2;
		fh3 = 0.25*ft3;
		fh4 = 0.25*ft4;
		fh5 = 0.25*ft5;

		F_x_wall[0][i][j][k]=fh1;
		F_x_wall[1][i][j][k]=fh2;
		F_x_wall[2][i][j][k]=fh3;
		F_x_wall[3][i][j][k]=fh4;
		if(Mtrl_kind!=3){
		F_x_wall[4][i][j][k]=fh5;
		}
		else{
		F_x_wall[4][i][j][k] = 0.0;
		}
	}
	else
	{
	//obtain the eigen vectors
	Real eigen_l[Emax][Emax], eigen_r[Emax][Emax], eigen_value[Emax];
	RoeAverage_x(i, j, k, eigen_l, eigen_r, eigen_value);
	//define type of artificial viscosity
	Real artificial_viscosity[Emax];
	for(int n=0; n<Emax; n++)
	{
		Real eigen_local_max = 0.0;
		for(int m=0; m<stencil_size; m++)
			eigen_local_max = AMAX1(eigen_local_max,  fabs(eigen_local[n][m+i-stencil_P][j][k]));
		artificial_viscosity[n] = Roe_type*eigen_value[n] + LLF_type*eigen_local_max + GLF_type*eigen_block[n];
	}
	//Scalar u, f, f_plus, f_minus
	//construct the right value & the left value scalar equations by characteristic reduction
	// at i+1/2 in x direction
	for(int n=0; n<Emax; n++){
		//characteristic decomposition
		for(int m=0; m<stencil_size; m++){
			uf[m] = 0.0;
			ff[m] = 0.0;
			for(int n1=0; n1<Emax; n1++){
				uf[m] += UI[n1][m+i-stencil_P][j][k]*eigen_l[n][n1];
				ff[m] += F_x[n1][m+i-stencil_P][j][k]*eigen_l[n][n1];
			}
			pp[m] = 0.5*(ff[m] + artificial_viscosity[n]*uf[m]);
			mm[m] = 0.5*(ff[m] - artificial_viscosity[n]*uf[m]);
		}
		// calculate the scalar numerical flux at x direction
		if(Mtrl_kind == 1 && CAVITATION==1 && (rho[i][j][k]<0.99995775*rho0 || rho[i+1][j][k]<0.99995775*rho0)){
			f_flux	= upwind_P(&pp[stencil_P], dx)
					+ upwind_M(&mm[stencil_P], dx);
		}
		else{
			f_flux = (*CurrentInterpolation_P)(&pp[stencil_P], dx)
				   + (*CurrentInterpolation_M)(&mm[stencil_P], dx);
		}

		// get Fp
		for(int n1=0; n1<Emax; n1++)
			_p[n][n1] = f_flux*eigen_r[n1][n];
	}
	// reconstruction the F_x terms
	for(int n=0; n<Emax; n++){
		F_x_wall[n][i][j][k] = 0.0;
		for(int n1=0; n1<Emax; n1++)
			F_x_wall[n][i][j][k] += _p[n1][n];
	}
}
}

int indicator(Real *f)
{
	using namespace std;
	int k;
	Real v0, v7,v1, v2, v3, v4, v5, v6;
	Real b1, b2, b3, b4, b5, b6, b7;
	Real a1, a2, a3, a4, a5, a6, a7;

	//assign value to v1, v2,...
	k = 0;
	v0 = *(f + k - 3);
	v1 = *(f + k - 2);
	v2 = *(f + k - 1);
	v3 = *(f + k);
	v4 = *(f + k + 1);
	v5 = *(f + k + 2);
	v6 = *(f + k + 3);
	v7 = *(f + k + 4);

	Real s1 = 13.0 / 12.0*(v2 - 2.0*v3 + v4)*(v2 - 2.0*v3 + v4) + 3.0 / 12.0*(v2 - v4)*(v2 - v4);
	Real s2 = 13.0 / 12.0*(v3 - 2.0*v4 + v5)*(v3 - 2.0*v4 + v5) + 3.0 / 12.0*(3.0*v3 - 4.0*v4 + v5)*(3.0*v3 - 4.0*v4 + v5);
	Real s3 = 13.0 / 12.0*(v1 - 2.0*v2 + v3)*(v1 - 2.0*v2 + v3) + 3.0 / 12.0*(v1 - 4.0*v2 + 3.0*v3)*(v1 - 4.0*v2 + 3.0*v3);
	Real s4 = fabs(10.0*v6*v6-49.0*v5*v6+29.0*v4*v6+61.0*v5*v5-73.0*v4*v5+22.0*v4*v4)/3.;  
	Real s5 = fabs(22.0*v2*v2-73.0*v1*v2+29.0*v0*v2+61.0*v1*v1-49.0*v0*v1+10.0*v0*v0)/3.;
	Real s6 = fabs(22.0*v7*v7-103.0*v6*v7+59.0*v5*v7+121.0*v6*v6-139.0*v5*v6+40.0*v5*v5)/3.;
		
	Real eps_hybrid;
	if (DIM_X + DIM_Y + DIM_Z == 1)
		eps_hybrid = 1.0e-4;
	else
	 	eps_hybrid = 1.0e-2;

	a1 = 1.0 / powern(s1 + eps_hybrid, 6.0);
	a2 = 1.0 / powern(s2 + eps_hybrid, 6.0);
	a3 = 1.0 / powern(s3 + eps_hybrid, 6.0);
	a4 = 1.0 / powern(s4 + eps_hybrid, 6.0);
	a5 = 1.0 / powern(s5 + eps_hybrid, 6.0);
	a6 = 1.0 / powern(s6 + eps_hybrid, 6.0);

	b1 = a1 / (a1 + a2 + a3 + a4+ a5 + a6);
	b2 = a2 / (a1 + a2 + a3 + a4+ a5 + a6);
	b3 = a3 / (a1 + a2 + a3 + a4+ a5 + a6);
	b4 = a4 / (a1 + a2 + a3 + a4+ a5 + a6);
	b5 = a5 / (a1 + a2 + a3 + a4+ a5 + a6);
	b6 = a6 / (a1 + a2 + a3 + a4+ a5 + a6);


	Real C_T = 1.e-6;

	int indicator_x;
	if (b1 > C_T && b2 > C_T && b3 > C_T && b4 > C_T&& b5 > C_T && b6 > C_T)
		indicator_x = 0;
	else
		indicator_x = 1;

	return indicator_x;
}


double CurrentInterpolation_P(double *f, double delta)
{
	using namespace std;
	int k;
	double v0, v1, v2, v3, v4, v5, v6, v7;
	double b1, b2, b3, b4, b5, b6;
	double a1, a2, a3, a4, a5, a6;
	double w1, w2, w3, w4, w5, w6;
	double d1, d2, d3, d4, d5, d6;
	double var1, var2, var3, var4, var5, var6, var7;
	double Variation1, Variation2, Variation3, Variation4, Variation5, Variation6;

	//assign value to v1, v2,...
	k = 0;
	v0 = *(f + k - 3);
	v1 = *(f + k - 2);
	v2 = *(f + k - 1);
	v3 = *(f + k);
	v4 = *(f + k + 1);
	v5 = *(f + k + 2);
	v6 = *(f + k + 3);
	v7 = *(f + k + 4);
	//smoothness indicator
	double s1 = 13.0 / 12.0*(v2 - 2.0*v3 + v4)*(v2 - 2.0*v3 + v4) + 3.0 / 12.0*(v2 - v4)*(v2 - v4);
	double s2 = 13.0 / 12.0*(v3 - 2.0*v4 + v5)*(v3 - 2.0*v4 + v5) + 3.0 / 12.0*(3.0*v3 - 4.0*v4 + v5)*(3.0*v3 - 4.0*v4 + v5);
	double s3 = 13.0 / 12.0*(v1 - 2.0*v2 + v3)*(v1 - 2.0*v2 + v3) + 3.0 / 12.0*(v1 - 4.0*v2 + 3.0*v3)*(v1 - 4.0*v2 + 3.0*v3);
	double s4 = 1. / 240.*fabs((2107.0*v3*v3 - 9402.0*v3*v4 + 11003.0*v4*v4 + 7042.0*v3*v5 - 17246.0*v4*v5 + 7043.0*v5*v5 - 1854.0*v3*v6 + 4642.0*v4*v6 - 3882.0*v5*v6 + 547.0*v6*v6));
	double s5 = 1. / 240.*fabs(v0*(547.0*v0 - 3882.0*v1 + 4642.0*v2 - 1854.0*v3) + v1*(7043.0*v1 - 17246.0*v2 + 7042.0*v3) + v2*(11003.0*v2 - 9402.0*v3) + 2107.0*v3*v3);
	double s6 = 1. / 5040.*fabs(v3*(107918.0*v3 - 649501.0*v4 + 758823.0*v5 - 411487.0*v6 + 86329.0*v7)
		+ v4*(1020563.0*v4 - 2462076.0*v5 + 1358458.0*v6 - 288007.0*v7)
		+ v5*(1521393.0*v5 - 1704396.0*v6 + 364863.0*v7)
		+ v6*(482963.0*v6 - 208501.0*v7)
		+ 22658.0*v7*v7);

	double tau8 = 1. / 62270208000.0 *fabs(v7*(75349098471.0*v7 - 1078504915264.0*v6 + 3263178215782.0*v5 - 5401061230160.0*v4 + 5274436892970.0*v3 - 3038037798592.0*v2 + 956371298594.0 *v1 - 127080660272.0*v0)
		+ v6*(3944861897609.0*v6 - 24347015748304.0*v5 + 41008808432890.0*v4 - 40666174667520.0*v3 + 23740865961334.0*v2 - 7563868580208.0 *v1 + 1016165721854.0*v0)
		+ v5*(38329064547231.0*v5 - 131672853704480.0*v4 + 132979856899250.0*v3 - 78915800051952.0*v2 + 25505661974314.0 *v1 - 3471156679072.0*v0)
		+ v4*(115451981835025.0*v4 - 238079153652400.0*v3 + 144094750348910.0*v2 - 47407534412640.0 *v1 + 6553080547830.0*v0)
		+ v3*(125494539510175.0*v3 - 155373333547520.0*v2 + 52241614797670.0 *v1 - 7366325742800.0*v0)
		+ v2*(49287325751121.0*v2 - 33999931981264.0 *v1 + 4916835566842.0*v0)
		+ v1*(6033767706599.0 *v1 - 1799848509664.0*v0)
		+ 139164877641.0*v0*v0);


	a1 = powern(1.0 + tau8 / (s1 + 1.0e-40), 6.0);
	a2 = powern(1.0 + tau8 / (s2 + 1.0e-40), 6.0);
	a3 = powern(1.0 + tau8 / (s3 + 1.0e-40), 6.0);
	a4 = powern(1.0 + tau8 / (s4 + 1.0e-40), 6.0);
	a5 = powern(1.0 + tau8 / (s5 + 1.0e-40), 6.0);
	a6 = powern(1.0 + tau8 / (s6 + 1.0e-40), 6.0);

	double inverse = 1. / (a1 + a2 + a3 + a4 + a5 + a6);

	b1 = a1 * inverse;
	b2 = a2 * inverse;
	b3 = a3 * inverse;
	b4 = a4 * inverse;
	b5 = a5 * inverse;
	b6 = a6 * inverse;
	// adaptive CT
	var2 = fabs(v2 - v1);
	var3 = fabs(v3 - v2);
	var4 = fabs(v4 - v3);
	var5 = fabs(v5 - v4);
	var6 = fabs(v6 - v5);

	double r_c = 0.23;
	double kethe = 1.e-3;
	double eps = 0.9*r_c*kethe*kethe / (1. - 0.9*r_c);

	double r_1 = (2.*var2*var3 + eps) / (var2*var2 + var3*var3 + eps);
	double r_2 = (2.*var3*var4 + eps) / (var3*var3 + var4*var4 + eps);
	double r_3 = (2.*var4*var5 + eps) / (var4*var4 + var5*var5 + eps);
	double r_4 = (2.*var5*var6 + eps) / (var5*var5 + var6*var6 + eps);

	double r_min = AMIN1(AMIN1(r_1, r_2), AMIN1(r_3, r_4));
	double delta_x = 1. - AMIN1(r_min / r_c, 1.);

	double incom_int = 10.5;
	double diff_int = incom_int - 7.;
	double decay = powern((1.0 - delta_x), 4.0) * (1.0 + 4.0 * delta_x);
	double power_int = incom_int - diff_int * (1. - decay);
	double cut_off = 1. / powern(10., power_int);
	//sharp selection
	b1 = b1 < cut_off ? 0. : 1.;
	b2 = b2 < cut_off ? 0. : 1.;
	b3 = b3 < cut_off ? 0. : 1.;
	b4 = b4 < cut_off ? 0. : 1.;
	b5 = b5 < cut_off ? 0. : 1.;
	b6 = b6 < cut_off ? 0. : 1.;

	Variation1 = -1.0 / 6.0*v2 + 5.0 / 6.0*v3 + 2.0 / 6.0*v4 - v3;
	Variation2 = 2. / 6.*v3 + 5. / 6.*v4 - 1. / 6.*v5 - v3;
	Variation3 = 2. / 6.*v1 - 7. / 6.*v2 + 11. / 6.*v3 - v3;
	Variation4 = 3. / 12.*v3 + 13. / 12.*v4 - 5. / 12.*v5 + 1. / 12.*v6 - v3;
	Variation5 = -3. / 12.*v0 + 13. / 12.*v1 - 23. / 12.*v2 + 25. / 12.*v3 - v3;
	Variation6 = 12. / 60.*v3 + 77. / 60.*v4 - 43. / 60.*v5 + 17. / 60.*v6 - 3. / 60.*v7 - v3;

	a1 = .4130855804023061  * b1;
	a2 = .2193140179474727  * b2;
	a3 = .06459052081827904 * b3;
	a4 = .1439236310125986  * b4;
	a5 = .02746675592227204 * b5;
	a6 = .1316194938970745  * b6;

	inverse = 1. / (a1 + a2 + a3 + a4 + a5 + a6);

	w1 = a1 * inverse;
	w2 = a2 * inverse;
	w3 = a3 * inverse;
	w4 = a4 * inverse;
	w5 = a5 * inverse;
	w6 = a6 * inverse;

	return	v3 + w1*Variation1
		+ w2*Variation2
		+ w3*Variation3
		+ w4*Variation4
		+ w5*Variation5
		+ w6*Variation6;
}


double CurrentInterpolation_M(double *f, double delta)
{
	using namespace std;
	int k;
	double v0, v1, v2, v3, v4, v5, v6, v7;
	double b1, b2, b3, b4, b5, b6;
	double a1, a2, a3, a4, a5, a6;
	double w1, w2, w3, w4, w5, w6;
	double d1, d2, d3, d4, d5, d6;
	double var1, var2, var3, var4, var5, var6, var7;
	double Variation1, Variation2, Variation3, Variation4, Variation5, Variation6;

	k = 1;
	v0 = *(f + k + 3);
	v1 = *(f + k + 2);
	v2 = *(f + k + 1);
	v3 = *(f + k);
	v4 = *(f + k - 1);
	v5 = *(f + k - 2);
	v6 = *(f + k - 3);
	v7 = *(f + k - 4);

	double s1 = 13.0 / 12.0*(v2 - 2.0*v3 + v4)*(v2 - 2.0*v3 + v4) + 3.0 / 12.0*(v2 - v4)*(v2 - v4);
	double s2 = 13.0 / 12.0*(v3 - 2.0*v4 + v5)*(v3 - 2.0*v4 + v5) + 3.0 / 12.0*(3.0*v3 - 4.0*v4 + v5)*(3.0*v3 - 4.0*v4 + v5);
	double s3 = 13.0 / 12.0*(v1 - 2.0*v2 + v3)*(v1 - 2.0*v2 + v3) + 3.0 / 12.0*(v1 - 4.0*v2 + 3.0*v3)*(v1 - 4.0*v2 + 3.0*v3);
	double s4 = 1. / 240.*fabs((2107.0*v3*v3 - 9402.0*v3*v4 + 11003.0*v4*v4 + 7042.0*v3*v5 - 17246.0*v4*v5 + 7043.0*v5*v5 - 1854.0*v3*v6 + 4642.0*v4*v6 - 3882.0*v5*v6 + 547.0*v6*v6));
	double s5 = 1. / 240.*fabs(v0*(547.0*v0 - 3882.0*v1 + 4642.0*v2 - 1854.0*v3) + v1*(7043.0*v1 - 17246.0*v2 + 7042.0*v3) + v2*(11003.0*v2 - 9402.0*v3) + 2107.0*v3*v3);
	double s6 = 1. / 5040.*fabs(v3*(107918.0*v3 - 649501.0*v4 + 758823.0*v5 - 411487.0*v6 + 86329.0*v7)
		+ v4*(1020563.0*v4 - 2462076.0*v5 + 1358458.0*v6 - 288007.0*v7)
		+ v5*(1521393.0*v5 - 1704396.0*v6 + 364863.0*v7)
		+ v6*(482963.0*v6 - 208501.0*v7)
		+ 22658.0*v7*v7);

	double tau8 = 1. / 62270208000.0 *fabs(v7*(75349098471.0*v7 - 1078504915264.0*v6 + 3263178215782.0*v5 - 5401061230160.0*v4 + 5274436892970.0*v3 - 3038037798592.0*v2 + 956371298594.0 *v1 - 127080660272.0*v0)
		+ v6*(3944861897609.0*v6 - 24347015748304.0*v5 + 41008808432890.0*v4 - 40666174667520.0*v3 + 23740865961334.0*v2 - 7563868580208.0 *v1 + 1016165721854.0*v0)
		+ v5*(38329064547231.0*v5 - 131672853704480.0*v4 + 132979856899250.0*v3 - 78915800051952.0*v2 + 25505661974314.0 *v1 - 3471156679072.0*v0)
		+ v4*(115451981835025.0*v4 - 238079153652400.0*v3 + 144094750348910.0*v2 - 47407534412640.0 *v1 + 6553080547830.0*v0)
		+ v3*(125494539510175.0*v3 - 155373333547520.0*v2 + 52241614797670.0 *v1 - 7366325742800.0*v0)
		+ v2*(49287325751121.0*v2 - 33999931981264.0 *v1 + 4916835566842.0*v0)
		+ v1*(6033767706599.0 *v1 - 1799848509664.0*v0)
		+ 139164877641.0*v0*v0);


	a1 = powern(1.0 + tau8 / (s1 + 1.0e-40), 6.0);
	a2 = powern(1.0 + tau8 / (s2 + 1.0e-40), 6.0);
	a3 = powern(1.0 + tau8 / (s3 + 1.0e-40), 6.0);
	a4 = powern(1.0 + tau8 / (s4 + 1.0e-40), 6.0);
	a5 = powern(1.0 + tau8 / (s5 + 1.0e-40), 6.0);
	a6 = powern(1.0 + tau8 / (s6 + 1.0e-40), 6.0);

	double inverse = 1. / (a1 + a2 + a3 + a4 + a5 + a6);

	b1 = a1 * inverse;
	b2 = a2 * inverse;
	b3 = a3 * inverse;
	b4 = a4 * inverse;
	b5 = a5 * inverse;
	b6 = a6 * inverse;

	var2 = fabs(v2 - v1);
	var3 = fabs(v3 - v2);
	var4 = fabs(v4 - v3);
	var5 = fabs(v5 - v4);
	var6 = fabs(v6 - v5);

	double r_c = 0.23;
	double kethe = 1.e-3;
	double eps = 0.9*r_c*kethe*kethe / (1. - 0.9*r_c);

	double r_1 = (2.*var2*var3 + eps) / (var2*var2 + var3*var3 + eps);
	double r_2 = (2.*var3*var4 + eps) / (var3*var3 + var4*var4 + eps);
	double r_3 = (2.*var4*var5 + eps) / (var4*var4 + var5*var5 + eps);
	double r_4 = (2.*var5*var6 + eps) / (var5*var5 + var6*var6 + eps);

	double r_min = AMIN1(AMIN1(r_1, r_2), AMIN1(r_3, r_4));
	double delta_x = 1. - AMIN1(r_min / r_c, 1.);

	double incom_int = 10.5;
	double diff_int = incom_int - 7.;
	double decay = powern((1.0 - delta_x), 4.0) * (1.0 + 4.0 * delta_x);
	double power_int = incom_int - diff_int * (1. - decay);
	double cut_off = 1. / powern(10., power_int);

	b1 = b1 < cut_off ? 0. : 1.;
	b2 = b2 < cut_off ? 0. : 1.;
	b3 = b3 < cut_off ? 0. : 1.;
	b4 = b4 < cut_off ? 0. : 1.;
	b5 = b5 < cut_off ? 0. : 1.;
	b6 = b6 < cut_off ? 0. : 1.;

	Variation1 = -1.0 / 6.0*v2 + 5.0 / 6.0*v3 + 2.0 / 6.0*v4 - v3;
	Variation2 = 2. / 6.*v3 + 5. / 6.*v4 - 1. / 6.*v5 - v3;
	Variation3 = 2. / 6.*v1 - 7. / 6.*v2 + 11. / 6.*v3 - v3;
	Variation4 = 3. / 12.*v3 + 13. / 12.*v4 - 5. / 12.*v5 + 1. / 12.*v6 - v3;
	Variation5 = -3. / 12.*v0 + 13. / 12.*v1 - 23. / 12.*v2 + 25. / 12.*v3 - v3;
	Variation6 = 12. / 60.*v3 + 77. / 60.*v4 - 43. / 60.*v5 + 17. / 60.*v6 - 3. / 60.*v7 - v3;

	a1 = .4130855804023061  * b1;
	a2 = .2193140179474727  * b2;
	a3 = .06459052081827904 * b3;
	a4 = .1439236310125986  * b4;
	a5 = .02746675592227204 * b5;
	a6 = .1316194938970745  * b6;

	inverse = 1. / (a1 + a2 + a3 + a4 + a5 + a6);

	w1 = a1 * inverse;
	w2 = a2 * inverse;
	w3 = a3 * inverse;
	w4 = a4 * inverse;
	w5 = a5 * inverse;
	w6 = a6 * inverse;

	return	v3 + w1*Variation1
		+ w2*Variation2
		+ w3*Variation3
		+ w4*Variation4
		+ w5*Variation5
		+ w6*Variation6;
}