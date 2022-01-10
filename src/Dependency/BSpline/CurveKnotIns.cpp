#include "CurveKnotIns.h"



CurveKnotIns::CurveKnotIns()
{
}


CurveKnotIns::~CurveKnotIns()
{
}

// Insert a knot r times.
CurveKnotIns::CurveKnotIns(int &np, int &p, std::vector<double> &UP, std::vector<Point4> &Pw, double &u, int &k, int &s, int &r, int &nq, std::vector<double> &UQ, std::vector<Point4> &Qw)
{
	/*Compute new curve from knot insertion  */
	/* Input: np, p, UP, Pw, u, k, s, r	*/
		/* np: number of contorl points	(-1)		  */
		/* p : degree								  */
		/* UP: Kont vector							  */
		/* Pw: contorl points in homogeneous form	  */
		/* u : new kont								  */
		/* k : u in [u_k, u_{k+1})					  */
		/* s : initial multiplicity of u			  */
		/* r : insert u r times						  */
	/* Output : nq, UQ, Qw */
		/* nq: number of contorl points				  */
		/* UQ: new kont vector						  */
		/* Qw: new contorl points in homogeneous form */
	nq = np + r;

	/* Load new knot vector */
	UQ.clear();
	UQ.resize(UP.size() + r);
	int mp = np + p + 1;
	for (int i = 0; i <= k; i++) UQ[i] = UP[i];
	for (int i = 1; i <= r; i++) UQ[k + i] = u;
	for (int i = k + 1; i <= mp; i++) UQ[i + r] = UP[i];

	/* Save unaltered control points */
	Qw.clear();
	Qw.resize(nq + 1);
	std::vector<Point4> Rw(p - s + 1);
	for (int i = 0; i <= k - p; i++) Qw[i] = Pw[i];
	for (int i = k - s; i <= np; i++) Qw[i + r] = Pw[i];
	for (int i = 0; i <= p - s; i++) Rw[i] = Pw[k - p + i];

	/* Insert the knot r times */
	int L;
	for (int j = 1; j <= r; j++)
	{
		L = k - p + j;
		for (int i = 0; i <= p - j - s; i++)
		{
			double alpha = (u - UP[L + i]) / (UP[i + k + 1] - UP[L + i]);
			Rw[i] = alpha * Rw[i + 1] + (1.0 - alpha) * Rw[i];
		}

		Qw[L] = Rw[0];
		Qw[k - s + r - j] = Rw[p - s - j];
	}
	for (int i = 1; i < p - s - r; i++)		/* Load remaining control points */
		Qw[L + i] = Rw[i];
}

// Insert many knots.
CurveKnotIns::CurveKnotIns(int &n, int &p, std::vector<double> &U, std::vector<Point4> &Pw, std::vector<double> &X, int &r, std::vector<double> &Ubar, std::vector<Point4> &Qw)
{
	/*Compute new curve from knot refinement */
	/* Input: n, p, U, Pw, X, r	***************** */
	/* n : number of contorl points	(-1)		  */
	/* p : degree								  */
	/* U : Kont vector							  */
	/* Pw: contorl points in homogeneous form	  */
	/* X : knot vector{x_0,...,x_r} to be inserted*/
	/* r : maximum index of new knot element x_i  */
	/* Output : Ubar, Qw ************************ */
	/* Ubar: new kont vector					  */
	/* Qw: new contorl points in homogeneous form */
	int a, b;
	int a_multi, b_multi;
	FindSpanMult(X[0], a, a_multi);
	FindSpanMult(X[r], b, b_multi);
	b = b + 1;

	/* initialize Ubar */
	Ubar.clear();
	Ubar.resize(U.size() + r + 1);
	for (int i = 0; i <= a; i++) Ubar[i] = U[i];
	for (int i = b; i <= n + p + 1; i++) Ubar[i + r + 1] = U[i];

	/* Save unaltered ctrl pts */
	Qw.clear();
	Qw.resize(Pw.size() + r + 1);
	for (int i = 0; i <= a - p; i++) Qw[i] = Pw[i];
	for (int i = b - 1; i <= n; i++) Qw[i + r + 1] = Pw[i];

	/* Insert the knots backward in X one by one */
	int i = b + p + 1, k = b + p + r;
	for (int j = r; j >= 0; j++)
	{
		while(X[j] <= U[i] && i > a)	// 为了得到初始的Qw(最多p+2个)，用于求插入 x_j∈(U_{i-1},U_i] 之后的ctlpts 
		{
			Qw[k - p - 1] = Pw[i - p - 1];
			Ubar[k] = U[i];
			k = k - 1;
			i = i - 1;
		}
		Qw[k - p - 1] = Pw[i - p - 1];
		for (int L = 1; L <= p; L++)
		{
			int index = k - p + 1;
			double beta = Ubar[k + L] - X[j]; // β = 1 - α
			if (abs(beta) < DBL_EPSILON) Qw[index - 1] = Qw[index];
			else
			{
				beta /= (Ubar[k + L] - U[i - p + 1]);
				Qw[index - 1] = beta * Qw[index - 1] + (1 - beta) * Qw[index];
			}		
		}
		Ubar[k] = X[j];
		k = k - 1;
	}
}

// Split Bspline to Bezier curves.
void CurveKnotIns::BsplineToBeziers_Curve(int &n, int &p, std::vector<double> U, std::vector<Point4> Pw, int &nb, std::vector<std::vector<Point4>> &Qw)
{
	/*Compute new curves from decomposition */
	/* Input: n, p, U, Pw *********************** */
	/* n : number of contorl points	(-1)		  */
	/* p : degree								  */
	/* U : Kont vector							  */
	/* Pw: contorl points in homogeneous form	  */
	/* Output : np, Qw ************************** */
	/* nb: number of Bezier sub-curve segment	  */
	/* Qw: new contorl points of Bezier sub-curves*/

	int m = n + p + 1;
	int a = p, b = p + 1;
	nb = 0;
	Qw.clear();
	for (int i = 0; i <= p; i++) Qw[nb][i] = Pw[i];

	while (b < m)
	{
		int multi;
		int i = b;
		while (b < m && ((U[b + 1] - U[b]) < DBL_EPSILON)) 
			b++;
		multi = b - i + 1;

		if (multi < p)
		{
			/* Compute alphas */		
			std::vector<double> alphas(p - multi);
			double numerator = U[b] - U[a];
			for (int j = 0; j < p - multi; j--)
				alphas[j] = numerator / (U[a + multi + 1 + j] - U[a]);

			/* Insert U[b] r times */
			int r = p - multi;
			for (int j = 1; j <= r; j++)
			{
				int save = r - j;
				int s = multi + j;
				for (int k = p; k >= s; k--)
				{
					double alpha = alphas[k - s];
					Qw[nb][k] = alpha * Qw[nb][k] + (1.0 - alpha) * Qw[nb][k - 1];
				}			
				if (b < m)
					Qw[nb + 1][save] = Qw[nb][p]; // control points of next segment
			}
		}
		nb++;

		/* Initialize for next segment */
		if (b < m)
		{
			for (int i = p - multi; i <= p; i++) 
				Qw[nb][i] = Pw[b - p + i];
			a = b;
			b = b + 1;
		}
	}
}


// Split Bezier to Bezier curves.
void CurveKnotIns::BezierDecomposition_Curve(int p, std::vector<Point4> &Pw, double u_min, double u_max, double u, std::vector<Point4> &left_Qw, std::vector<Point4> &right_Qw)
{
	/*Compute new curves from decomposition */
	/* Input: p, Pw, u_min, u_max, u ************ */
	/* p : degree								  */
	/* u_min : minimum parameter of current curve */
	/* u_max : maximum parameter of current curve */
	/* u : split the curve in u					  */
	/* Pw: contorl points in homogeneous form	  */
	/* Output : left_Qw, right_Qw *************** */
	/* left_Qw: contorl points of left segment	  */
	/* right_Qw: contorl points of right segment  */

	assert (u <= u_max && u >= u_min);

	left_Qw.clear();
	right_Qw.clear();
	left_Qw.resize(p + 1);
	right_Qw.resize(p + 1);
	std::vector<Point4> Rw = Pw;
	left_Qw[0] = Rw[0];
	right_Qw[p] = Rw[p];
	double alpha = (u - u_min) / (u_max - u_min);
	for (int j = 1; j <= p; j++)
	{
		for (int i = p; i >= j; i--)
		{
			Rw[i] = alpha * Rw[i] + (1.0 - alpha) * Rw[i - 1];
		}

		left_Qw[j] = Rw[j];
		right_Qw[p - j] = Rw[p];          /////////// 0515
	}
}


void CurveKnotIns::FindSpanMult(double u, int &knot_span, int &multi)
{
	auto knots = curve_.GetKnots();

	knot_span = 0;	
	int p = curve_.GetDegree();
	for (int i = p; i < knots.size() - p; i++)
	{
		if (u < knots[i])
		{
			knot_span = i - 1;
			break;
		}
	}

	multi = 0;
	for (int j = knot_span; j >= 0; j--)
	{
		if (abs(u - knots[j]) < DBL_EPSILON)
			multi++;
		else
			break;
	}
}