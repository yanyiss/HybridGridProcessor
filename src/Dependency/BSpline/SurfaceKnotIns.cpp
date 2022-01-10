#include "SurfaceKnotIns.h"



SurfaceKnotIns::SurfaceKnotIns()
{
}


SurfaceKnotIns::~SurfaceKnotIns()
{
}

// Insert a knot r times
SurfaceKnotIns::SurfaceKnotIns(int &np, int &p, std::vector<double> &UP, int &mp, int &q, std::vector<double> &VP, std::vector<std::vector<Point4>> &Pw, DIRECTION dir, int &uv, int &k, int &s, int &r, int &nq, std::vector<double> &UQ, int &mq, std::vector<double> &VQ, std::vector<std::vector<Point4>> &Qw)
{
	/* Surface knot insertion */
	/* Input: np, p, UP, mp, q, VP, Pw, dir, uv, k, s, r	*/
	/* np : number of contorl points in u direction(-1)	*/
	/* p  : degree in u direction						*/
	/* UP : Kont vector in u direction					*/
	/* mp : number of contorl points in v direction		*/
	/* q  : degree in v direction						*/
	/* VP : Kont vector in v direction					*/
	/* Pw : contorl points								*/
	/* dir: the knot insertion is implemented in U/V	*/
	/* uv : new kont									*/
	/* k  : uv in [u_k, u_{k+1}) or [v_k, v_{k+1})		*/
	/* s  : initial multiplicity of uv					*/
	/* r  : insert uv r times							*/
	/* Output : nq, UQ, mq, VQ, Qw */
	/* nq: number of new contorl points in u direction	*/
	/* UQ: new kont vector in u direction				*/
	/* mq: number of new contorl points in v direction	*/
	/* VQ: new kont vector in v direction				*/
	/* Qw: new contorl points							*/

	if (dir == U_DIRECTION)
	{
		nq = np + r;
		UQ.clear();
		UQ.resize(UP.size() + r);
		for (int i = 0; i <= k; i++) UQ[i] = UP[i];
		for (int i = 1; i <= r; i++) UQ[k + i] = uv;
		for (int i = k + 1; i <= np + p + 1; i++) UQ[i + r] = UP[i];

		mq = mp;
		VQ.clear();
		VQ = VP;

		/* Save the alphas */
		std::vector<std::vector<double>> alpha(p - s, std::vector<double>(r));
		for (int j = 1; j <= r; j++)
		{
			int L = k - p + j;
			for (int i = 0; i <= p - j - s; i++)
				alpha[i][j - 1] = (uv - UP[L + i]) / (UP[i + k + 1] - UP[i + L]);
		}

		Qw.clear();
		Qw.resize(nq + 1, std::vector<Point4>(mq + 1));
		for (int row = 0; row <= mp; row++)
		{
			/* Save unaltered control points */
			for (int i = 0; i <= k - p; i++) Qw[i][row] = Pw[i][row];
			for (int i = k - s; i <= np; i++) Qw[i + r][row] = Pw[i][row];

			/* Load auxiliary control points. */
			std::vector<Point4> Rw(p + 1 - s);
			for (int i = 0; i <= p - s; i++) Rw[i] = Pw[i + k - p][row];
			for (int j = 1; j <= r; j++)
			{
				for (int i = 0; i <= p - s - j; i++)
				{
					Rw[i] = alpha[i][j - 1] * Rw[i + 1] + (1 - alpha[i][j - 1]) * Rw[i];
				}
				Qw[k - p + j][row] = Rw[0];
				Qw[k - s + r - j][row] = Rw[p - s - j];
			}

			/* Load the remaining control points */
			for (int i = 1; i < p - s - r; i++)
			{
				Qw[k - p + r + i][row] = Rw[i];
			}
		}

	}

	else
	{
		nq = np;
		UQ.clear();
		UQ = UP;

		mq = mp + r;
		VQ.clear();
		VQ.resize(VP.size() + r);
		for (int j = 0; j <= k; j++) VQ[j] = VP[j];
		for (int j = 1; j <= r; j++) VQ[k + j] = uv;
		for (int j = k + 1; j <= mp + q + 1; j++) VQ[j + r] = VP[j];

		/* Save the alphas */
		std::vector<std::vector<double>> alpha(r, std::vector<double>(q - s));
		for (int i = 1; i <= r; i++)
		{
			int L = k - q + i;
			for (int j = 0; j <= q - i - s; j++)
				alpha[i - 1][j] = (uv - VP[L + j]) / (VP[j + k + 1] - VP[j + L]);
		}


		Qw.clear();
		Qw.resize(nq + 1, std::vector<Point4>(mq + 1));
		for (int col = 0; col <= np; col++)
		{
			/* Save unaltered control points */
			for (int j = 0; j <= k - q; j++) Qw[col][j] = Pw[col][j];
			for (int j = k - s; j <= mp; j++) Qw[col][j + r] = Pw[col][j];

			/* Load auxiliary control points. */
			std::vector<Point4> Rw(q + 1 - s);
			for (int j = 0; j <= q - s; j++) Rw[j] = Pw[col][j + k - q];
			for (int i = 1; i <= r; i++)
			{
				for (int j = 0; j <= q - s - i; j++)
				{
					Rw[j] = alpha[i - 1][j] * Rw[j + 1] + (1 - alpha[i - 1][j]) * Rw[j];
				}
				Qw[col][k - q + i] = Rw[0];
				Qw[col][k - s + r - i] = Rw[q - s - i];
			}

			/* Load the remaining control points */
			for (int j = 1; j < q - s - r; j++)
			{
				Qw[col][k - q + r + j] = Rw[j];
			}
		}
	}

}

// Insert many knots.
SurfaceKnotIns::SurfaceKnotIns(int &n, int &p, std::vector<double> &U, int &m, int &q, std::vector<double> &V, std::vector<std::vector<Point4>> &Pw, std::vector<double> &X, int &r, DIRECTION dir, std::vector<double> &Ubar, std::vector<double> &Vbar, std::vector<std::vector<Point4>> &Qw)
{
	/* Surface knot refinement */
	/* Input: n, p, U, m, q, V, Pw, X, r, dir ***********/
	/* n  : number of contorl points in u direction(-1)	*/
	/* p  : degree in u direction						*/
	/* U  : Kont vector in u direction					*/
	/* m  : number of contorl points in v direction		*/
	/* q  : degree in v direction						*/
	/* V  : Kont vector in v direction					*/
	/* Pw : contorl points								*/
	/* X  : knot vector{x_0,...,x_r} to be inserted		*/
	/* r  : insert uv r times							*/
	/* dir: the knot refinement is implemented in U/V	*/
	/* Output : Ubar, Vbar, Qw **************************/
	/* Ubar: new kont vector in u direction				*/
	/* Vbar: new kont vector in v direction				*/
	/* Qw: new contorl points							*/
	int a, b;
	int a_multi, b_multi;
	FindSpanMult(X[0], dir, a, a_multi);
	FindSpanMult(X[r], dir, b, b_multi);
	b = b + 1;

	if (dir == U_DIRECTION)
	{
		/* copy V into Vbar */
		Vbar.clear();
		Vbar.resize(V.size());
		Vbar = V;

		/* initialize Ubar */ 
		Ubar.clear();
		Ubar.resize(U.size() + r + 1);
		for (int i = 0; i <= a; i++) Ubar[i] = U[i];
		for (int i = b; i <= n + p + 1; i++) Ubar[i + r + 1] = U[i];

		/* Save unaltered ctrl pts */
		Qw.clear();
		Qw.resize(n + r + 2, std::vector<Point4>(m + 1));
		for (int row = 0; row <= m; row++)
		{
			for (int i = 0; i <= a - p; i++) Qw[i][row] = Pw[i][row];
			for (int i = b - 1; i <= n; i++) Qw[i + r + 1][row] = Pw[i][row];
		}
		
		/* Insert the knots backward in X one by one */
		int i = b + p + 1, k = b + p + r;
		for (int j = r; j >= 0; j--)
		{
			while (X[j] <= U[i] && i > a)
			{
				for (int row = 0; row <= m; row++) 
					Qw[k - p - 1][row] = Pw[i - p - 1] [row];
				k = k - 1;
				i = i - 1;
			}

			for (int row = 0; row <= m; row++)
				Qw[k - p - 1][row]= Qw[k - p][row];
			for (int L = 1; L <= p; L++)
			{
				int index = k - p + L;
				double beta = Ubar[k + L] - X[j]; // ¦Â = 1 - ¦Á
				if (abs(beta) < DBL_EPSILON)
				{
					for (int row = 0; row <= m; row++)
						Qw[index - 1][row] = Qw[index][row];
				}				
				else
				{
					beta /= (Ubar[k + L] - U[i - p + 1]);
					for (int row = 0; row <= m; row++)
						Qw[index - 1][row] = beta * Qw[index - 1][row] + (1.0 - beta) * Qw[index][row];
				}
			}
			Ubar[k] = X[j];
			k = k - 1;
		}		
	}
	
	else
	{
		/* Similar code as above with u and v directional parameters switched */

		/* copy U into Ubar */
		Ubar.clear();
		Ubar.resize(U.size());
		Ubar = U;

		/* initialize Vbar */
		Vbar.clear();
		Vbar.resize(V.size() + r + 1);
		for (int j = 0; j <= a; j++) Vbar[j] = V[j];
		for (int j = b; j <= m + q + 1; j++) Vbar[j + r + 1] = V[j];

		/* Save unaltered ctrlpts */
		Qw.clear();
		Qw.resize(n + 1, std::vector<Point4>(m + r + 2));
		for (int col = 0; col <= n; col++)
		{
			for (int j = 0; j <= a - q; j++) Qw[col][j] = Pw[col][j];
			for (int j = b - 1; j <= m; j++) Qw[col][j + r + 1] = Pw[col][j];
		}

		/* Insert the knots backward in X one by one */
		int j = b + q + 1, k = b + q + r;
		for (int i = r; i >= 0; i--)
		{
			while (X[i] <= V[j] && j > a)
			{
				for (int col = 0; col <= n; col++)
					Qw[col][k - q - 1] = Pw[col][j - q - 1];
				k = k - 1;
				j = j - 1;
			}

			for (int col = 0; col <= n; col++)
				Qw[col][k - q - 1] = Qw[col][k - q];

			for (int L = 1; L <= q; L++)
			{
				int index = k - q + L;
				double beta = Vbar[k + L] - X[i]; // ¦Â = 1 - ¦Á
				if (abs(beta) < DBL_EPSILON)
				{
					for (int col = 0; col <= n; col++)
						Qw[col][index - 1] = Qw[col][index];
				}
				else
				{
					beta /= (Vbar[k + L] - V[j - q + 1]);
					for (int col = 0; col <= n; col++)
						Qw[col][index - 1] = beta * Qw[col][index - 1] + (1.0 - beta) * Qw[col][index];
				}
			}
			Vbar[k] = X[i];
			k = k - 1;
		}
	}

}

// Split Bspline to Bezier surfaces.
void SurfaceKnotIns::BsplineToBeziers_Surface(int &n, int &p, std::vector<double> &U, int &m, int &q, std::vector<double> &V, std::vector<std::vector<Point4>> &Pw, DIRECTION &dir, int &nb, std::vector<std::vector<std::vector<Point4>>>Qw)
{
	/* Surface knot refinement */
	/* Input: n, p, U, m, q, V, Pw, dir *****************/
	/* n  : number of contorl points in u direction(-1)	*/
	/* p  : degree in u direction						*/
	/* U  : Kont vector in u direction					*/
	/* m  : number of contorl points in v direction		*/
	/* q  : degree in v direction						*/
	/* V  : Kont vector in v direction					*/
	/* Pw : contorl points								*/
	/* dir: the decomposition is implemented in U/V		*/
	/* Output : nb, Qw **********************************/
	/* nb: number of Bezier segment in one direction	*/
	/* Qw: new contorl points							*/

	Qw.clear();
	nb = 0;

	if (dir == U_DIRECTION)
	{		
		int a = p, b = p + 1;

		for (int i = 0; i <= p; i++)
			for (int row = 0; row <= m; row++)
				Qw[nb][i][row] = Pw[i][row];

		while (b < n + p + 1)
		{
			int multi;
			int i = b;
			while (b < n + p + 1 && ((U[b + 1] - U[b]) < DBL_EPSILON))
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
						for (int row = 0; row <= m; row++)
							Qw[nb][k][row] = alpha * Qw[nb][k][row] + (1.0 - alpha) * Qw[nb][k - 1][row];
					}

					if (b < n + p + 1)
						for (int row = 0; row <= m; row++)
							Qw[nb + 1][save][row] = Qw[nb][p][row]; // control points of next segment
				}
			}
			nb++;

			if (b < n + p + 1)
			{
				for (int i = p - multi; i <= p; i++)
					for (int row = 0; row <= m; row++)
						Qw[nb][i][row] = Pw[b - p + i][row];
					
				a = b;
				b = b + 1;
			}

		}
	}

	else
	{
		int a = q, b = q + 1;

		for (int j = 0; j <= q; j++)
			for (int col = 0; col <= n; col++)
				Qw[nb][col][j] = Pw[col][j];

		while (b < m + q + 1)
		{
			int multi;
			int i = b;
			while (b < m + q + 1 && ((V[b + 1] - V[b]) < DBL_EPSILON))
				b++;
			multi = b - i + 1;

			if (multi < q)
			{
				/* Compute alphas */
				std::vector<double> alphas(q - multi);
				double numerator = V[b] - V[a];
				for (int i = 0; i < q - multi; i--)
					alphas[i] = numerator / (V[a + multi + 1 + i] - V[a]);

				/* Insert U[b] r times */
				int r = q - multi;
				for (int i = 1; i <= r; i++)
				{
					int save = r - i;
					int s = multi + i;
					for (int k = q; k >= s; k--)
					{
						double alpha = alphas[k - s];
						for (int col = 0; col <= n; col++)
							Qw[nb][col][k] = alpha * Qw[nb][col][k] + (1.0 - alpha) * Qw[nb][col][k - 1];
					}

					if (b < m + q + 1)
						for (int col = 0; col <= n; col++)
							Qw[nb + 1][col][save] = Qw[nb][col][p]; // control points of next segment
				}
			}
			nb++;

			if (b < m + q + 1)
			{
				for (int j = q - multi; j <= q; j++)
					for (int col = 0; col <= n; col++)
						Qw[nb][col][j] = Pw[col][b - q + j];

				a = b;
				b = b + 1;
			}

		}
	}

}

// Split Bezier to Bezier suefaces.
void SurfaceKnotIns::BezierDecomposition_Surface(int &p, int &q, std::vector<std::vector<Point4>> Pw, double u_min, double u_max, double v_max, double v_min, double uv, DIRECTION &dir, std::vector<std::vector<Point4>> &left_Qw, std::vector<std::vector<Point4>> &right_Qw)
{
	/* Compute new surfaces from decomposition */
	/* Input: p, q, Pw, u_min, u_max, v_min, v_max, dir */
	/* p  : degree in u direction						*/
	/* q  : degree in v direction						*/
	/* Pw : contorl points								*/
	/* u_min : minimum parameter in u direction			*/
	/* u_max : maximum parameter in u direction			*/
	/* v_min : minimum parameter in v direction			*/
	/* v_max : maximum parameter in v direction			*/
	/* u : split the surface in uv						*/
	/* dir: the decomposition is implemented in U/V		*/
	/* Output : left_Qw, right_Qw ********************* */
	/* left_Qw: contorl points of left segment	  */
	/* right_Qw: contorl points of right segment  */

	left_Qw.clear();
	right_Qw.clear();
	left_Qw.resize(p + 1, std::vector<Point4>(q + 1));
	right_Qw.resize(p + 1, std::vector<Point4>(q + 1));
	std::vector<std::vector<Point4>> Rw = Pw;

	if (dir == U_DIRECTION)
	{
		if (uv <= u_min || uv >= u_max) return;

//		left_Qw[0] = Rw[0];
		for (int row = 0; row <= q; row++)
		{
			left_Qw[0][row] = Rw[0][row];
			right_Qw[p][row] = Rw[p][row];

			double alpha = (uv - u_min) / (u_max - u_min);
			for (int j = 1; j <= p; j++)
			{
				for (int i = p; i >= j; i--)
				{
					Rw[i][row] = alpha * Rw[i][row] + (1.0 - alpha) * Rw[i - 1][row];
				}

				left_Qw[j][row] = Rw[j][row];
				right_Qw[p - j][row] = Rw[p][row];
			}
		}
	}

	else
	{
		if (uv <= v_min || uv >= v_max) return;

		for (int col = 0; col <= p; col++)
		{
			left_Qw[col][0] = Rw[col][0];
			right_Qw[col][q] = Rw[col][q];

			double alpha = (uv - v_min) / (v_max - v_min);
			for (int i = 1; i <= q; i++)
			{
				for (int j = q; j >= i; j--)
				{
					Rw[col][j] = alpha * Rw[col][j] + (1.0 - alpha) * Rw[col][j - 1];
				}

				left_Qw[col][i] = Rw[col][i];
				right_Qw[col][q - i] = Rw[col][q];
			}
		}
	}
}

void SurfaceKnotIns::FindSpanMult(double uv, DIRECTION dir, int &knot_span, int &multi)
{	
	int degree;
	std::vector<double>knots;
	if (dir == U_DIRECTION)
	{
		degree = surface_.GetDegreeU();
		knots = surface_.GetKnotsU();
	}					
	else
	{
		degree = surface_.GetDegreeV();
		knots = surface_.GetKnotsV();
	}
		
	knot_span = 0;
	for (int i = degree; i < knots.size()- degree; i++)
	{
		if (uv < knots[i])
		{
			knot_span = i - 1;
			break;
		}
	}

	multi = 0;
	for (int j = knot_span; j >= 0; j--)
	{
		if (abs(uv - knots[j]) < DBL_EPSILON)
			multi++;
		else
			break;
	}
}