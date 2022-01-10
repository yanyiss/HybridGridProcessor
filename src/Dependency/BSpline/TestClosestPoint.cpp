#include "TestClosestPoint.h"
#include <math.h>


TestClosestPoint::TestClosestPoint(BezierCurve &curve, Point &testpoint, double &dist, std::vector<ProjectionPointToCurve> &closest_pts)
{	
	P = testpoint;
	closest_pts.clear();
	closest_pts_c.clear();
	RemainingCurveSegments.clear();
	RemainingCurveSegments.push_back(curve);
	curr_min_dist = dist;
	
	while (!RemainingCurveSegments.empty())
	{
		auto curve_ = RemainingCurveSegments.front();

		bool isEndPoint = TestEndPointOfCurve(curve_, P);
		if(isEndPoint)  // the closet point is one of the end points
		{
			RemainingCurveSegments.pop_front();
			continue;
		}
		
		bool isUniqueInside = TestIntUniquenessOfCurve(curve_, P);
		if (isUniqueInside) // abadon or unique(conputed by Newton's)
		{
			RemainingCurveSegments.pop_front();
			continue;
		}
					
		// Subdivide to 2 sub-curves.
		BezierCurve left_curve, right_curve;
		double t_min, t_max;
		curve_.GetParameterSpan(t_min, t_max);
		curve_.Split((t_min + t_max) / 2.0, left_curve, right_curve);
		RemainingCurveSegments.push_back(left_curve);
		RemainingCurveSegments.push_back(right_curve);
		RemainingCurveSegments.pop_front();
	}

	dist = curr_min_dist;
	closest_pts = closest_pts_c;

	// �õ��� closest_pts_c ��Ϊ���
}

// ���һ�������߽���ĸ�����
// ȡ�����Ǹ������ڻ�����Ƿ�ö��ϵ������
// ������ǣ�����������
// ����ǣ��뵱ǰ����̾���Ƚϣ������������Ϣ
// ����������
bool TestClosestPoint::TestEndPointOfCurve(BezierCurve &curve, Point P)
{
	double t_min, t_max;
	double deg = curve.GetDegree();
	curve.GetParameterSpan(t_min, t_max);
	auto ctrlpts = curve.GetControlPoints();

	Point P_0 = ctrlpts[0], P_p = ctrlpts[deg];
	double dist_P0 = (P_0 - P).norm();
	double dist_Pp = (P_p - P).norm();
	UpdateCurveMessage(t_min, dist_P0);
	UpdateCurveMessage(t_max, dist_Pp);

	if (dist_P0 < dist_Pp)
	{
		Point VecPQ = P_0 - P;		
		for (int i = 1; i <= deg; i++)
		{
			if (VecPQ.dot(ctrlpts[i] - P_0) <= 0)
				return false;
		}	
	}
	else
	{
		Point VecPQ = P_p - P;
		for (int i = 0; i < deg; i++)
		{
			if (VecPQ.dot(ctrlpts[i] - P_p) <= 0)
				return false;
		}
	}

	return true;
}

// ����squared distance function
// ��ϵ��ȫ���ڵ�ǰ��̾��룬ֱ��������һ�Σ�����ϸ��
// ���򣬼��ϵ���Ƿ�ֻ��һ���ֲ���С
// ����Ψһ�Բ�����,�������أ�����ϸ��
// ���ǣ�ת��Newton����⣬����ϸ�֡�
bool TestClosestPoint::TestIntUniquenessOfCurve(BezierCurve &curve, Point P)
{
	double t_min, t_max;
	curve.GetParameterSpan(t_min, t_max);
	if (t_max - t_min > 1E-6)
	{
		SquaredDistanceFunc_Curve distFunc_curve(curve, P);
		std::vector<double> distFuncCoeffs = distFunc_curve.GetCoeffs();
		int deg_sdf = distFunc_curve.GetDegree();

		auto min_itr = std::min_element(distFuncCoeffs.begin(), distFuncCoeffs.end());
		if (*min_itr > pow(curr_min_dist, 2)) 		// the minimum coefficient is grater than alpha, abadon this segment
			return true;					// not unique, but don't need to subodivide this segment anymore

		int idx = std::distance(distFuncCoeffs.begin(), min_itr);
		if (idx == 0 || idx == deg_sdf) return false;

		for (int i = idx; i > 0; i--)
		{
			if (distFuncCoeffs[i] > distFuncCoeffs[i - 1])
				return false;
		}
		for (int i = idx; i < deg_sdf; i++)
		{
			if (distFuncCoeffs[i + 1] < distFuncCoeffs[i])
				return false;
		}
	}

	/*********************************************************/
	// ��Newton Raphson��������������㣬�뵱ǰ��̾���Ƚϣ�����
	//double t_min, t_max;
	//curve.GetParameterSpan(t_min, t_max);
	double optimal_t = (t_min + t_max) / 2.0;

	NewtonRaphson NR_Iteration(P, curve, optimal_t);

	double dist = (P - curve(optimal_t)).norm();
	UpdateCurveMessage(optimal_t, dist);
	/*********************************************************/

	return true;
}



TestClosestPoint::TestClosestPoint(BezierSurface &surface, Point &testpoint, std::vector<ProjectionPointToSurface> & closest_pts)
{
	P = testpoint;
	closest_pts.clear();
	closest_pts_s.clear();
	RemainingSurfacePatches.clear();
	RemainingSurfacePatches.push_back(surface);
	curr_min_dist = 1.0E10;

	while (!RemainingSurfacePatches.empty())
	{
		auto surface_ = RemainingSurfacePatches.front();

		bool isCornerPoint = TestCornerPointOfSurface(surface_, P);
		if (isCornerPoint)
		{
			RemainingSurfacePatches.pop_front();	 // the closet point is one of the coner points
			continue;
		}
		
		bool isOnBoundaryCurve = TestBoudaryCurveOfSurface(surface_, P);
		if (isOnBoundaryCurve)
		{
			RemainingSurfacePatches.pop_front();	 // the closet point is on boundary curve
			continue;
		}
		
		bool isUniqueInside = TestIntUniquenessOfSurface(surface_, P);
		if (isUniqueInside)
		{
			RemainingSurfacePatches.pop_front();// abadon or unique(conputed by Newton's)
			continue;
		}

		BezierSurface left_surface, right_surface;
		double u_min, u_max, v_min, v_max;
		surface_.GetParameterSpan(u_min, u_max, v_min, v_max);
		double u_len = u_max - u_min;
		double v_len = v_max - v_min;
		if (u_len > v_len)
		{
			surface_.Split((u_min + u_max) / 2.0, U_DIRECTION, left_surface, right_surface);
		}
		else
		{
			surface_.Split((v_min + v_max) / 2.0, V_DIRECTION, left_surface, right_surface);
		}
		RemainingSurfacePatches.push_back(left_surface);
		RemainingSurfacePatches.push_back(right_surface);
		RemainingSurfacePatches.pop_front();
	}

	closest_pts = closest_pts_s;
	// �õ��� closest_pts_s ��Ϊ���	
}

// ���һ���ĸ��ǵ��ĸ�����
// ȡ�����Ǹ������ڻ�����Ƿ��patch�ϵ������
// ������ǣ�����������
// ����ǣ��뵱ǰ����̾���Ƚϣ�
// �����ڵ�ǰ��С���룬����¼����С�ڣ���ɾ�������µ�ǰ�������Ϣ��
// ����������
bool TestClosestPoint::TestCornerPointOfSurface(BezierSurface &surface, Point P)
{
	auto u_deg = surface.GetDegreeU();
	auto v_deg = surface.GetDegreeV();
	auto ctrlpts = surface.GetControlPoints();
	double u_min, u_max, v_min, v_max;
	surface.GetParameterSpan(u_min, u_max, v_min, v_max);
	
	Point P_00 = ctrlpts[0][0], P_0q = ctrlpts[0][v_deg];
	Point P_p0 = ctrlpts[u_deg][0], P_pq = ctrlpts[u_deg][v_deg];
	double d_00 = (P_00 - P).norm();
	double d_0q = (P_0q - P).norm();
	double d_p0 = (P_p0 - P).norm();
	double d_pq = (P_pq - P).norm();
	UpdateSurfaceMessage(u_min, v_min, d_00);
	UpdateSurfaceMessage(u_min, v_max, d_0q);
	UpdateSurfaceMessage(u_max, v_min, d_p0);
	UpdateSurfaceMessage(u_max, v_max, d_pq);

	// P_00 is the closet one of 4 corner points
	if (d_00 <= d_0q && d_00 <= d_p0 && d_00 <= d_pq)
	{
		Point VecPQ = P_00 - P;
		for (int j = 1; j <= v_deg; j++)
		{
			if (VecPQ.dot(ctrlpts[0][j] - P_00) < 0)
				return false;
		}
		for (int i = 1; i <= u_deg; i++)
		{
			for (int j = 0; j <= v_deg; j++)
			{
				if (VecPQ.dot(ctrlpts[i][j] - P_00) < 0)
					return false;
			}
		}
	}

	// P_0q is the closet one of 4 corner points
	else if(d_0q <= d_00 && d_0q <= d_p0 && d_0q <= d_pq)
	{
		Point VecPQ = P_0q - P;
		for (int j = 0; j < v_deg; j++)
		{
			if (VecPQ.dot(ctrlpts[0][j] - P_0q) < 0)
				return false;
		}
		for (int i = 1; i <= u_deg; i++)
		{
			for (int j = 0; j <= v_deg; j++)
			{
				if (VecPQ.dot(ctrlpts[i][j] - P_0q) < 0)
					return false;
			}
		}
	}

	// P_p0 is the closet one of 4 corner points
	else if (d_p0 <= d_00 && d_p0 <= d_0q && d_p0 <= d_pq)
	{
		Point VecPQ = P_p0 - P;
		for (int j = 1; j <= v_deg; j++)
		{
			if (VecPQ.dot(ctrlpts[u_deg][j] - P_p0) < 0)
				return false;
		}
		for (int i = 0; i < u_deg; i++)
		{
			for (int j = 0; j <= v_deg; j++)
			{
				if (VecPQ.dot(ctrlpts[i][j] - P_p0) < 0)
					return false;
			}
		}
	}

	// P_pq is the closet one of 4 corner points
	else if (d_pq <= d_00 && d_pq <= d_0q && d_pq <= d_p0)
	{
		Point VecPQ = P_pq - P;
		for (int j = 0; j < v_deg; j++)
		{
			if (VecPQ.dot(ctrlpts[u_deg][j] - P_pq) < 0)
				return false;
		}
		for (int i = 0; i < u_deg; i++)
		{
			for (int j = 0; j <= v_deg; j++)
			{
				if (VecPQ.dot(ctrlpts[i][j] - P_pq) < 0)
					return false;
			}
		}
	}
	
	return true;
}

// �ÿ��ƶ���difference���ڻ��ж�������Ƿ���4���߽�������
// 4���߽�: u = u_0, u_p; v = v_0, v_q �ֱ���һ��
// ��ĳ���߽��Ϸ��������ͶӰ�����������ת��Ϊ���ڸ������ڲ��������ͶӰ����
// ��4���߽������Ͼ����������������ڸ�patch�ϵ������������ڲ�
// ת����֤�ڲ���Ψһ��
bool TestClosestPoint::TestBoudaryCurveOfSurface(BezierSurface &surface, Point P)
{
	auto u_deg = surface.GetDegreeU();
	auto v_deg = surface.GetDegreeV();
	auto ctrlpts = surface.GetControlPoints();
	double u_min, u_max, v_min, v_max;
	surface.GetParameterSpan(u_min, u_max, v_min, v_max);

	bool isOnBdrCurve = true;

	// S(u_0,v)
	for (int j = 0; j <= v_deg && isOnBdrCurve; j++)
	{
		auto VecPQ = ctrlpts[0][j] - P;
		for (int i = 0; i < u_deg && isOnBdrCurve; i++)
		{
			for (int l = 0; l <= v_deg; l++)
			{
				if (VecPQ.dot(ctrlpts[i + 1][l] - ctrlpts[i][l]) < 0)
				{
					isOnBdrCurve = false;
					break;				
				}
			}
		}
	}
	
	if (isOnBdrCurve == true)
	{
		/**************����߽��ϵ������**************/
		std::vector<Point> bdr_ctrlpts(ctrlpts[0]);	
		BezierCurve curve(v_deg, v_min, v_max, bdr_ctrlpts);
		bool isRational = surface.GetRational();
		if (isRational)
		{
			curve.SetRational(true);
			auto w = surface.GetWeights();
			curve.SetWeights(w[0]);
		}
		
		double min_dist = curr_min_dist;
		TestClosestPoint(curve, P, curr_min_dist, closest_pts_c);
		if (curr_min_dist < min_dist)
			closest_pts_s.clear();

		for (int i = 0; i < closest_pts_c.size(); i++)
		{
			auto proj = closest_pts_c[i];
			UpdateSurfaceMessage(u_min, proj.t, proj.dist);
		}
	
		return true;
	}

	// S(u_p, v)
	isOnBdrCurve = true;
	for (int j = 0; j <= v_deg && isOnBdrCurve; j++)
	{
		auto VecPQ = P - ctrlpts[u_deg][j];
		for (int i = 0; i < u_deg && isOnBdrCurve; i++)
		{
			for (int l = 0; l <= v_deg; l++)
			{
				if (VecPQ.dot(ctrlpts[i + 1][l] - ctrlpts[i][l]) < 0)
				{
					isOnBdrCurve = false;
					break;
				}				
			}
		}
	}
	
	if (isOnBdrCurve == true) // the closet point is on boudary curve S(u_p,v)
	{
		/**************����߽��ϵ������**************/
		std::vector<Point> bdr_ctrlpts(ctrlpts[u_deg]);
		BezierCurve curve(v_deg, v_min, v_max, bdr_ctrlpts);
		bool isRational = surface.GetRational();
		if (isRational)
		{
			curve.SetRational(true);
			auto w = surface.GetWeights();
			curve.SetWeights(w[u_deg]);
		}
	
		double min_dist = curr_min_dist;
		TestClosestPoint(curve, P, curr_min_dist, closest_pts_c);
		if (curr_min_dist < min_dist)
			closest_pts_s.clear();

		for (int i = 0; i < closest_pts_c.size(); i++)
		{
			auto proj = closest_pts_c[i];
			UpdateSurfaceMessage(u_max, proj.t, proj.dist);
		}
	
		return true;
	}

	// S(u,v_0)
	isOnBdrCurve = true;
	for (int i = 0; i <= u_deg && isOnBdrCurve; i++)
	{
		auto VecPQ = ctrlpts[i][0] - P;
		for (int k = 0; k <= u_deg && isOnBdrCurve; k++)
		{
			for (int j = 0; j < v_deg; j++)
			{
				if (VecPQ.dot(ctrlpts[k][j + 1] - ctrlpts[k][j]) < 0)
				{
					isOnBdrCurve = false;
					break;
				}
			}
		}
	}
	
	if (isOnBdrCurve == true) // the closet point is on boudary curve S(u,v_0)
	{
		/**************����߽��ϵ������**************/
		std::vector<Point> bdr_ctrlpts(u_deg + 1);
		for (int i = 0; i <= u_deg; i++)
		{
			bdr_ctrlpts[i] = ctrlpts[i][0];
		}
	
		BezierCurve curve(u_deg, u_min, u_max, bdr_ctrlpts);
		bool isRational = surface.GetRational();
		if (isRational)
		{
			curve.SetRational(true);
			auto w = surface.GetWeights();
			std::vector<double> weights(u_deg + 1);
			for (int i = 0; i < u_deg; i++)
			{
				weights[i] = w[i][0];
			}
			curve.SetWeights(weights);
		}
		
		double min_dist = curr_min_dist;
		TestClosestPoint(curve, P, curr_min_dist, closest_pts_c);
		if (curr_min_dist < min_dist)
			closest_pts_s.clear();

		for (int i = 0; i < closest_pts_c.size(); i++)
		{
			auto proj = closest_pts_c[i];
			UpdateSurfaceMessage(proj.t, v_min, proj.dist);
		}
	
		return true;
	}

	// S(u, v_q)
	isOnBdrCurve = true;
	for (int i = 0; i <= u_deg && isOnBdrCurve; i++)
	{
		auto VecPQ = P - ctrlpts[i][v_deg];
		for (int k = 0; k <= u_deg; k++)
		{
			for (int j = 0; j < v_deg && isOnBdrCurve; j++)
			{
				if (VecPQ.dot(ctrlpts[k][j + 1] - ctrlpts[k][j]) < 0)
				{
					isOnBdrCurve = false;
					break;
				}
			}
		}
	}
	if (isOnBdrCurve == true) // the closet point is on boudary curve S(u,v_q)
	{
		/**************����߽��ϵ������**************/
		std::vector<Point> bdr_ctrlpts(u_deg + 1);
		for (int i = 0; i <= u_deg; i++)
		{
			bdr_ctrlpts[i] = ctrlpts[i][v_deg];
		}
	
		BezierCurve curve(u_deg, u_min, u_max, bdr_ctrlpts);
		bool isRational = surface.GetRational();
		if (isRational)
		{
			curve.SetRational(true);
			auto w = surface.GetWeights();
			std::vector<double> weights(u_deg + 1);
			for (int i = 0; i < u_deg; i++)
			{
				weights[i] = w[i][v_deg];
			}
			curve.SetWeights(weights);
		}

		double min_dist = curr_min_dist;
		TestClosestPoint(curve, P, curr_min_dist, closest_pts_c);
		if (curr_min_dist < min_dist)
			closest_pts_s.clear();

		for (int i = 0; i < closest_pts_c.size(); i++)
		{
			auto proj = closest_pts_c[i];
			UpdateSurfaceMessage(proj.t, v_max, proj.dist);
		}

		return true;
	}

	return false;
}

bool TestClosestPoint::TestIntUniquenessOfSurface(BezierSurface &surface, Point P, bool)
{
	double u_min, u_max, v_min, v_max;
	surface.GetParameterSpan(u_min, u_max, v_min, v_max);

	// ����һ��ƽ��(�����ǵ�)�� �������п��ƶ��㵽ƽ��������룬С��һ����ֵʱ������ӽ���ƽ��
	// ��ʱ����Newton�����
	auto ctrlpts = surface.GetControlPoints();
	auto P_00 = ctrlpts[0][0];
	auto P_p0 = ctrlpts.back()[0];
	auto P_0q = ctrlpts[0].back();
	auto normal = (P_p0 - P_00).cross(P_0q - P_00);
	double ave_dist = ((P_p0 - P_00).norm() + (P_0q - P_00).norm()) / 2;
	double max_dis = ave_dist;
	
	for (int i = 0; i < ctrlpts.size(); i++)
	{
		for (int j = 0; j < ctrlpts[0].size(); j++)
		{
			double dist = normal.dot(ctrlpts[i][j]) / normal.norm();
			max_dis = std::min(max_dis, dist);
		}
	}

	if (max_dis < 0.05 * ave_dist || std::max(u_max -u_min, v_max - v_min) < 1E-6)
	{
		double u_min, u_max, v_min, v_max;
		surface.GetParameterSpan(u_min, u_max, v_min, v_max);
		double optimal_u = (u_min + u_max) / 2.0;
		double optimal_v = (v_min + v_max) / 2.0;

		NewtonRaphson NR_Iteration(P, surface, optimal_u, optimal_v);

		double dist = (P - surface(optimal_u, optimal_v)).norm();
		UpdateSurfaceMessage(optimal_u, optimal_v, dist);

		return true;
	}
	else
	{
		return false;
	}
}

bool TestClosestPoint::TestIntUniquenessOfSurface(BezierSurface &surface, Point P)
{
	double u_min, u_max, v_min, v_max;
	surface.GetParameterSpan(u_min, u_max, v_min, v_max);
	if (u_max - u_min < 1E-6 || v_max - v_min < 1E-6)
	{
		double u_min, u_max, v_min, v_max;
		surface.GetParameterSpan(u_min, u_max, v_min, v_max);
		double optimal_u = (u_min + u_max) / 2.0;
		double optimal_v = (v_min + v_max) / 2.0;

		NewtonRaphson NR_Iteration(P, surface, optimal_u, optimal_v);

		double dist = (P - surface(optimal_u, optimal_v)).norm();
		UpdateSurfaceMessage(optimal_u, optimal_v, dist);

		return true;
	}

	if (surface.GetRational())
		return false;
	/*************** orthogonal projection: F = G = 0(neccessary)******************/
	/*************** compute the ggradient coefficients of F & G ******************/
	int p = surface.GetDegreeU(); // u_degree
	int q = surface.GetDegreeV(); // v_degree
	auto ctrlpts = surface.GetControlPoints();
	std::vector<double> binominal_p = Binomial(p);
	std::vector<double> binominal_p_ = Binomial(p - 1);
	std::vector<double> binominal_p__ = Binomial(p - 2);
	std::vector<double> binominal_q = Binomial(q);
	std::vector<double> binominal_q_ = Binomial(q - 1);
	std::vector<double> binominal_q__ = Binomial(q - 2);
	std::vector<double> binominal_2p = Binomial(2 * p);
	std::vector<double> binominal_2p_ = Binomial(2 * p - 1);
	std::vector<double> binominal_2p__ = Binomial(2 * p - 2);
	std::vector<double> binominal_2q = Binomial(2 * q);
	std::vector<double> binominal_2q_ = Binomial(2 * q - 1);
	std::vector<double> binominal_2q__ = Binomial(2 * q - 2);

	std::vector<std::vector<double>> F1_coeffs(2 * p - 1, std::vector<double>(2 * q + 1));
	std::vector<std::vector<double>> F2_coeffs(2 * p, std::vector<double>(2 * q));
	std::vector<std::vector<double>> G1_coeffs(2 * p, std::vector<double>(2 * q));
	std::vector<std::vector<double>> G2_coeffs(2 * p + 1, std::vector<double>(2 * q - 1));

	for (int k = 0; k <= 2 * p - 2; k++)
	{
		for (int l = 0; l <= 2 * q; l++)
		{
			double temp_coeff = 0.0;
			for (int i = std::max(0, k - p + 1); i <= std::min(k, p - 1); i++)
			{
				for (int j = std::max(0, l - q); j <= std::min(l, q); j++)
				{
					temp_coeff += (ctrlpts[i + 1][j].dot(ctrlpts[k - i + 1][l - j]) - 2 * ctrlpts[i + 1][j].dot(ctrlpts[k - i][l - j]) + ctrlpts[i][j].dot(ctrlpts[k - i][l - j])) * binominal_p_[i] * binominal_p_[k - i] * binominal_q[j] * binominal_q[l - j];
				}
			}
			F1_coeffs[k][l] = p * p * temp_coeff / binominal_2p__[k] / binominal_2q[l];

			temp_coeff = 0.0;
			for (int i = std::max(0, k - p); i <= std::min(k, p - 2); i++)
			{
				for (int j = std::max(0, l - q); j <= std::min(l, q); j++)
				{
					temp_coeff += (ctrlpts[i + 2][j] - 2 * ctrlpts[i + 1][j] + ctrlpts[i][j]).dot(ctrlpts[k - i][l - j] - P) * binominal_p__[i] * binominal_p[k - i] * binominal_q[j] * binominal_q[l - j];
				}
			}
			F1_coeffs[k][l] += p * (p - 1) * temp_coeff / binominal_2p__[k] / binominal_2q[l];

		}
	}

	for (int k = 0; k <= 2 * p - 1; k++)
	{
		for (int l = 0; l <= 2 * q - 1; l++)
		{
			double temp_coeff = 0.0;
			for (int i = std::max(0, k - p); i <= std::min(k, p - 1); i++)
			{
				for (int j = std::max(0, l - q + 1); j <= std::min(l, q); j++)
				{
					temp_coeff += (ctrlpts[i + 1][j].dot(ctrlpts[k - i][l - j + 1]) - ctrlpts[i + 1][j].dot(ctrlpts[k - i][l - j]) - ctrlpts[i][j].dot(ctrlpts[k - i][l - j + 1]) + ctrlpts[i][j].dot(ctrlpts[k - i][l - j])) * binominal_p_[i] * binominal_p[k - i] * binominal_q[j] * binominal_q_[l - j];
				}
			}
			
			for (int i = std::max(0, k - p); i <= std::min(k, p - 1); i++)
			{
				for (int j = std::max(0, l - q); j <= std::min(l, q - 1); j++)
				{
					temp_coeff += (ctrlpts[i + 1][j + 1] - ctrlpts[i + 1][j] - ctrlpts[i][j + 1] + ctrlpts[i][j]).dot(ctrlpts[k - i][l - j] - P) * binominal_p_[i] * binominal_p[k - i] * binominal_q_[j] * binominal_q[l - j];
				}
			}
			F2_coeffs[k][l] = p * q * temp_coeff / binominal_2p_[k] / binominal_2q_[l];

		}
	}

	G1_coeffs = F2_coeffs;

	for (int k = 0; k <= 2 * p; k++)
	{
		for (int l = 0; l <= 2 * q - 2; l++)
		{
			double temp_coeff = 0.0;
			for (int i = std::max(0, k - p); i <= std::min(k, p); i++)
			{
				for (int j = std::max(0, l - q + 1); j <= std::min(l, q - 1); j++)
				{
					temp_coeff += (ctrlpts[i][j + 1].dot(ctrlpts[k - i][l - j + 1]) - 2 * ctrlpts[i][j + 1].dot(ctrlpts[k - i][l - j]) + ctrlpts[i][j].dot(ctrlpts[k - i][l - j])) * binominal_p[i] * binominal_p[k - i] * binominal_q_[j] * binominal_q_[l - j];
				}
			}
			G2_coeffs[k][l] = q * q * temp_coeff / binominal_2p[k] / binominal_2q__[l];

			temp_coeff = 0.0;
			for (int i = std::max(0, k - p); i <= std::min(k, p); i++)
			{
				for (int j = std::max(0, l - q); j <= std::min(l, q - 2); j++)
				{
					temp_coeff += (ctrlpts[i][j + 2] - 2 * ctrlpts[i][j + 1] + ctrlpts[i][j]).dot(ctrlpts[k - i][l - j] - P) * binominal_p[i] * binominal_p[k - i] * binominal_q__[j] * binominal_q[l - j];
				}
			}
			G2_coeffs[k][l] += q * (q - 1) * temp_coeff / binominal_2p[k] / binominal_2q__[l];

		}
	}

	/**************************** degree devation ****************************/
	std::vector<std::vector<double>> f1(2 * p, std::vector<double>(2 * q + 1));
	std::vector<std::vector<double>> f2(2 * p, std::vector<double>(2 * q + 1));
	std::vector<std::vector<double>> g1(2 * p + 1, std::vector<double>(2 * q));
	std::vector<std::vector<double>> g2(2 * p + 1, std::vector<double>(2 * q));

	for (int l = 0; l <= 2 * q; l++)
	{
		f1[0][l] = F1_coeffs[0][l];
		f1[2 * p - 1] = F1_coeffs[2 * p - 2];
		for (int k = 1; k < 2 * p - 1; k++)
		{
			f1[k][l] = ((2 * p - 1 - k) * F1_coeffs[k][l] + k * F1_coeffs[k - 1][l]) / (2 * p - 1);
		}
	}

	for (int k = 0; k <= 2 * p - 1; k++)
	{
		f2[k][0] = F2_coeffs[k][0];
		f2[k][2 * q] = F2_coeffs[k][2 * q - 1];
		for (int l = 1; l < 2 * q; l++)
		{
			f2[k][l] = ((2 * q - l) * F2_coeffs[k][l] + l * F2_coeffs[k][l - 1]) / (2 * q);
		}
	}

	for (int l = 0; l <= 2 * q - 1; l++)
	{
		g1[0][l] = G1_coeffs[0][l];
		g1[2 * p] = G1_coeffs[2 * p - 1];
		for (int k = 1; k < 2 * p; k++)
		{
			g1[k][l] = ((2 * p - k) * G1_coeffs[k][l] + k * G1_coeffs[k - 1][l]) / (2 * p);
		}
	}

	for (int k = 0; k <= 2 * p; k++)
	{
		g2[k][0] = G2_coeffs[k][0];
		g2[k][2 * q - 1] = G2_coeffs[k][2 * q - 2];
		for (int l = 1; l < 2 * q - 1; l++)
		{
			g2[k][l] = ((2 * q - 1 - l) * G2_coeffs[k][l] + l * G2_coeffs[k][l - 1]) / (2 * q - 1);
		}
	}

	/***************************** normal cone ******************************/
	// axis 
	// ���⣺Ҫ��Ҫ��ÿ��ϵ�������ȹ�һ����ȡƽ��
	std::vector<Eigen::Vector2d> gradF(2 * p * (2 * q + 1));
	std::vector<Eigen::Vector2d> gradG(2 * q * (2 * p + 1));
	for (int i = 0; i < 2 * p; i++)
	{
		for (int j = 0; j <= 2 * q; j++)
		{
			gradF.push_back(Eigen::Vector2d(f1[i][j], f2[i][j]));
		}
	}

	for (int i = 0; i <= 2 * p; i++)
	{
		for (int j = 0; j < 2 * q; j++)
		{
			gradG.push_back(Eigen::Vector2d(g1[i][j], g2[i][j]));
		}
	}

	Eigen::Vector2d axis_F(0.0, 0.0);
	Eigen::Vector2d axis_G(0.0, 0.0);
	for (auto itr = gradF.begin(); itr != gradF.end(); itr++)
	{
		axis_F += *itr;
	}
	for (auto itr = gradG.begin(); itr != gradG.end(); itr++)
	{
		axis_G += *itr;
	}
	axis_F /= gradF.size();
	axis_G /= gradG.size();

	double norm_f = axis_F.norm();
	double norm_g = axis_G.norm();
	double angle = acos(axis_F.dot(axis_G) / norm_f / norm_g);
	double alpha = 0.0, beta = 0.0;
	for (auto itr = gradF.begin(); itr != gradF.end(); itr++)
	{
		double normal = (*itr).norm();
		alpha = std::max(alpha, acos(axis_F.dot(*itr) / norm_f / normal));
	}
	for (auto itr = gradG.begin(); itr != gradG.end(); itr++)
	{
		double normal = (*itr).norm();
		beta = std::max(beta, acos(axis_F.dot(*itr) / norm_g / normal));
	}

	if (angle < M_PI - alpha - beta)
		return false;
	else
	{
		double u_min, u_max, v_min, v_max;
		surface.GetParameterSpan(u_min, u_max, v_min, v_max);
		double optimal_u = (u_min + u_max) / 2.0;
		double optimal_v = (v_min + v_max) / 2.0;

		NewtonRaphson NR_Iteration(P, surface, optimal_u, optimal_v);

		double dist = (P - surface(optimal_u, optimal_v)).norm();
		UpdateSurfaceMessage(optimal_u, optimal_v, dist);

		return true;
	}
}


TestClosestPoint::~TestClosestPoint()
{

}

void TestClosestPoint::UpdateCurveMessage(double t, double dist)
{
	if (abs(dist - curr_min_dist) < 1E-8)
	{
		for (auto &tmp: closest_pts_c)
		{
			if (abs(t - tmp.t) < 1E-4)	//�ظ�
				return;
		}
		closest_pts_c.push_back(ProjectionPointToCurve{ t, dist });
	}
	else 
	{
		if (dist < curr_min_dist) 
		{
			closest_pts_c.clear();
			curr_min_dist = dist;
			closest_pts_c.push_back(ProjectionPointToCurve{ t, dist });
		}		
	}
}

void TestClosestPoint::UpdateSurfaceMessage(double u, double v, double dist)
{
	if (abs(dist - curr_min_dist) < 1E-8)
	{
		for (auto &tmp : closest_pts_s)
		{
			if (abs(u - tmp.u) < 1E-4 && abs(v - tmp.v) < 1E-4)	//�ظ�
				return;
		}
		closest_pts_s.push_back(ProjectionPointToSurface{ u, v, dist });
	}
	else 
	{
		if (dist < curr_min_dist)
		{
			closest_pts_s.clear();
			curr_min_dist = dist;
			closest_pts_s.push_back(ProjectionPointToSurface{ u, v, dist });
		}		
	}
}


TestClosestPoint::TestClosestPoint(BSplineCurve &curve, Point &testpoint, double &dist, std::vector<ProjectionPointToCurve> & closest_pts)
{
	P = testpoint;
	closest_pts.clear();
	closest_pts_c.clear();
	RemainingCurveSegments.clear();
	std::vector<BezierCurve> bezier_curves;
	curve.BSplineToBezier(bezier_curves);
	for (auto &c : bezier_curves)
	{
		RemainingCurveSegments.push_back(c);
	}
	curr_min_dist = dist;

	while (!RemainingCurveSegments.empty())
	{
		auto curve_ = RemainingCurveSegments.front();

		bool isEndPoint = TestEndPointOfCurve(curve_, P);
		if (isEndPoint)  // the closet point is one of the end points
		{
			RemainingCurveSegments.pop_front();
			continue;
		}

		bool isUniqueInside = TestIntUniquenessOfCurve(curve_, P);
		if (isUniqueInside) // abadon or unique(conputed by Newton's)
		{
			RemainingCurveSegments.pop_front();
			continue;
		}

		// Subdivide to 2 sub-curves.
		BezierCurve left_curve, right_curve;
		double t_min, t_max;
		curve_.GetParameterSpan(t_min, t_max);
		curve_.Split((t_min + t_max) / 2.0, left_curve, right_curve);
		RemainingCurveSegments.push_back(left_curve);
		RemainingCurveSegments.push_back(right_curve);
		RemainingCurveSegments.pop_front();
	}

	dist = curr_min_dist;
	closest_pts = closest_pts_c;
}

TestClosestPoint::TestClosestPoint(BSplineSurface &surface, Point &testpoint, std::vector<ProjectionPointToSurface> & closest_pts)
{
	P = testpoint;
	closest_pts.clear();
	closest_pts_s.clear();
	RemainingSurfacePatches.clear();
	std::vector<BezierSurface> bezier_surfaces;
	surface.BSplineToBezier(bezier_surfaces);
	for (auto &s : bezier_surfaces)
	{
		RemainingSurfacePatches.push_back(s);
	}	
	curr_min_dist = 1.0E10;

	while (!RemainingSurfacePatches.empty())
	{
		auto surface_ = RemainingSurfacePatches.front();

		bool isCornerPoint = TestCornerPointOfSurface(surface_, P);
		if (isCornerPoint)
		{
			RemainingSurfacePatches.pop_front();	 // the closet point is one of the coner points
			continue;
		}

		bool isOnBoundaryCurve = TestBoudaryCurveOfSurface(surface_, P);
		if (isOnBoundaryCurve)
		{
			RemainingSurfacePatches.pop_front();	 // the closet point is on boundary curve
			continue;
		}
		
		bool isUniqueInside = TestIntUniquenessOfSurface(surface_, P, true);
		if (isUniqueInside)
		{
			RemainingSurfacePatches.pop_front();// abadon or unique(conputed by Newton's)
			continue;
		}

		BezierSurface left_surface, right_surface;
		double u_min, u_max, v_min, v_max;
		surface_.GetParameterSpan(u_min, u_max, v_min, v_max);
		double u_len = u_max - u_min;
		double v_len = v_max - v_min;
		if (u_len > v_len)
		{
			surface_.Split((u_min + u_max) / 2.0, U_DIRECTION, left_surface, right_surface);
		}
		else
		{
			surface_.Split((v_min + v_max) / 2.0, V_DIRECTION, left_surface, right_surface);
		}
		RemainingSurfacePatches.push_back(left_surface);
		RemainingSurfacePatches.push_back(right_surface);
		RemainingSurfacePatches.pop_front();
	}

	closest_pts = closest_pts_s;
}


// ��û�п��������ͶӰ��ֹ������һ���߽�������

// ϸ��֮������ܶ࣬��������һ��hyperplane clipping / KDOP
// ����������ʱ

// ���ߡ�����Ĵ���Ϊ1ʱ����2�׵����������