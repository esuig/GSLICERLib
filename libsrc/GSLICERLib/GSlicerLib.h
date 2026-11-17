
#pragma once

#include <sstream>

namespace gslicerlib
{
	namespace constants
	{
		const double MPI_HALF = 1.5707963267948966192313216916398;
		const double MPI = 3.1415926535897932384626433832795;
		const double MPI_2 = 6.283185307179586476925286766559;
		const double GET_TO_RAD = (MPI / 180);
		const double RAD_TO_DEG = (180 / MPI);

		const double COINC_THRESH = 1.0e-10;
	}


	/**
	* @brief Represents a three-dimensional point or vector.
	*
	* CPoint3 provides basic 3D vector operations, including arithmetic,
	* length calculations, normalization, dot product, cross product,
	* and distance calculations.
	*/
	class CPoint3
	{
	public:

		/**
		* @brief Union providing access to the vector coordinates.
		*
		* Allows the three coordinates to be accessed either through
		* an array or through individually named members.
		*/
		union
		{
			/**
			* @brief Array access to the vector coordinates (x, y, z).
			*/
			double m_pc[3];

			struct
			{
				/** @brief X coordinate */
				double m_x;

				/** @brief Y coordinate */
				double m_y;

				/** @brief Z coordinate */
				double m_z;
			};
		};

		/**
		* @brief Default constructor.
		*
		* Initializes all coordinates to zero.
		*/
		CPoint3()
		{
			m_x = 0.0;
			m_y = 0.0;
			m_z = 0.0;
		}

		/**
		* @brief Constructs a point from three coordinates.
		*
		* @param f_x Initial X coordinate.
		* @param f_y Initial Y coordinate.
		* @param f_z Initial Z coordinate.
		*/
		CPoint3(const double& f_x,
				const double& f_y,
				const double& f_z)
		{
			m_pc[0] = f_x;
			m_pc[1] = f_y;
			m_pc[2] = f_z;
		}

		/**
		* @brief Copy constructor.
		*
		* Initializes this point using the coordinates of another point.
		*
		* @param f_p Point to copy.
		*/
		CPoint3(const CPoint3& f_p)
		{
			m_pc[0] = f_p.m_pc[0];
			m_pc[1] = f_p.m_pc[1];
			m_pc[2] = f_p.m_pc[2];
		}

		/**
		* @brief Calculates the squared length of the vector.
		*
		* Computes the sum of the squared coordinates without performing
		* a square root operation.
		*
		* @return Squared vector length.
		*/
		inline double LengthSq() const
		{
			return m_x * m_x + m_y * m_y + m_z * m_z;
		}

		/**
		* @brief Calculates the length (magnitude) of the vector.
		*
		* @return Euclidean vector length.
		*/
		inline double Length() const
		{
			return std::sqrt(LengthSq());
		}

		/**
		* @brief Normalizes this vector in place.
		*
		* Divides each coordinate by the vector length, resulting in
		* a unit vector with a length of one.
		*
		* If the vector has zero length, its coordinates remain unchanged.
		*
		* @return Reference to this vector after normalization.
		*/
		inline CPoint3& UnitizeInPlace()
		{
			const double l_length = Length();
			if (l_length > 0.0)
			{
				m_x /= l_length;
				m_y /= l_length;
				m_z /= l_length;
			}

			return *this;
		}

		/**
		* @brief Returns a normalized copy of this vector.
		*
		* Creates a unit vector without modifying the original instance.
		* If the vector has zero length, an unchanged copy is returned.
		*
		* @return Normalized copy of this vector.
		*/
		inline CPoint3 GetUnitVector() const
		{
			const double l_length = Length();
			if (l_length > 0.0)
			{
				return CPoint3(m_x / l_length, m_y / l_length, m_z / l_length);
			}

			return *this;
		}

		/**
		* @brief Adds another vector to this vector.
		*
		* @param f_p Vector to add.
		* @return Resulting vector.
		*/
		inline CPoint3 operator+(const CPoint3& f_p) const
		{
			return CPoint3(m_x + f_p.m_x, m_y + f_p.m_y, m_z + f_p.m_z);
		}

		/**
		* @brief Subtracts another vector from this vector.
		*
		* @param f_p Vector to subtract.
		* @return Resulting vector.
		*/
		inline CPoint3 operator-(const CPoint3& f_p) const
		{
			return CPoint3(m_x - f_p.m_x, m_y - f_p.m_y, m_z - f_p.m_z);
		}

		/**
		* @brief Multiplies this vector by a scalar.
		*
		* @param f_scalar Scalar multiplication factor.
		* @return Resulting scaled vector.
		*/
		inline CPoint3 operator*(double f_scalar) const
		{
			return CPoint3(m_x * f_scalar, m_y * f_scalar, m_z * f_scalar);
		}

		/**
		* @brief Multiplies a vector by a scalar.
		*
		* Allows scalar multiplication with the scalar as the left operand.
		*
		* @param f_scalar Scalar multiplication factor.
		* @param f_p Vector to multiply.
		* @return Resulting scaled vector.
		*/
		inline friend CPoint3 operator*(double f_scalar, const CPoint3& f_p)
		{
			return CPoint3(f_p.m_x * f_scalar, f_p.m_y * f_scalar, f_p.m_z * f_scalar);
		}

	private:

	};

	/**
	* @brief Calculates the dot product of two vectors.
	*
	* @param f_a First vector.
	* @param f_b Second vector.
	* @return Dot product of the two vectors.
	*/
	inline static double Dot(const CPoint3& f_a, const CPoint3& f_b)
	{
		return f_a.m_x * f_b.m_x + f_a.m_y * f_b.m_y + f_a.m_z * f_b.m_z;
	}

	/**
	* @brief Calculates the cross product of two vectors.
	*
	* The resulting vector is perpendicular to both input vectors
	* and follows the right-hand rule.
	*
	* @param f_a First vector.
	* @param f_b Second vector.
	* @return Cross product of the two vectors.
	*/
	inline static CPoint3 Cross(const CPoint3& f_a, const CPoint3& f_b)
	{
		return CPoint3(f_a.m_y * f_b.m_z - f_a.m_z * f_b.m_y,
			f_a.m_z * f_b.m_x - f_a.m_x * f_b.m_z,
			f_a.m_x * f_b.m_y - f_a.m_y * f_b.m_x);
	}

	/**
	* @brief Calculates the squared distance between two points.
	*
	* Computes the squared Euclidean distance without performing
	* a square root operation.
	*
	* @param f_a First point.
	* @param f_b Second point.
	* @return Squared Euclidean distance between the two points.
	*/
	inline static double DistanceSq(const CPoint3& f_a, const CPoint3& f_b)
	{
		const double l_dx = f_b.m_x - f_a.m_x;
		const double l_dy = f_b.m_y - f_a.m_y;
		const double l_dz = f_b.m_z - f_a.m_z;

		return l_dx * l_dx + l_dy * l_dy + l_dz * l_dz;
	}

	/**
	* @brief Calculates the Euclidean distance between two points.
	*
	* @param f_a First point.
	* @param f_b Second point.
	* @return Euclidean distance between the two points.
	*/
	inline static double Distance(const CPoint3& f_a, const CPoint3& f_b)
	{
		return sqrt(DistanceSq(f_a, f_b));
	}


	/**
	* @class CIntersectionSegment
	* @brief Represents a line segment resulting from the intersection
	*        between a triangle and a plane
	*
	* Each segment is oriented according to the normal of the triangle
	* that generated it. The triangle normal points to the right-hand
	* side of the oriented segment.This orientation convention ensures 
	* that:
	* - External profiles are reconstructed counter-clockwise.
	* - Internal profiles (holes) are reconstructed clockwise.
	*/
	class CIntersectionSegment
	{
	public:

		/**
		* @brief Constructs an intersection segment with the specified endpoints
		*        and associated mesh information
		*
		* @param f_pt_start Starting point of the oriented intersection segment
		* @param f_pt_end Ending point of the oriented intersection segment
		* @param f_start_edge_id Identifier of the mesh edge associated with
		*                        the starting point. Set to -1 if the starting
		*                        point corresponds to a mesh vertex
		* @param f_end_edge_id Identifier of the mesh edge associated with
		*                      the ending point. Set to -1 if the ending
		*                      point corresponds to a mesh vertex
		* @param f_face_id Identifier of the mesh face (triangle) that generated
		*                  the intersection segment
		*/
		CIntersectionSegment(const CPoint3& f_pt_start,
							 const CPoint3& f_pt_end,
							 const uint64_t f_start_edge_id,
							 const uint64_t f_end_edge_id,
							 int f_face_id):
									m_intersection_start_pt(f_pt_start),
									m_intersection_end_pt(f_pt_end),
									m_start_edge_id(f_start_edge_id),
									m_end_edge_id(f_end_edge_id),
									m_face_id(f_face_id),
									m_segment_is_used(false)
		{

		}

		/**
		* @brief Default constructor.
		*
		* Initializes the segment endpoints to the lowest representable
		* double value in each coordinate, the edge identifiers to zero,
		* the face identifier to -1, and the usage flag to false.
		*
		* The initialized endpoints can be used to identify an uninitialized
		* intersection segment.
		*/
		CIntersectionSegment():
					m_intersection_start_pt(std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest()),
					m_intersection_end_pt(std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest()),
					m_start_edge_id(0),
					m_end_edge_id(0),
					m_face_id(-1),
					m_segment_is_used(false)
		{

		}

		/**
		* @brief Returns the starting point of the intersection segment.
		*
		* @return Constant reference to the starting point.
		*/
		inline const CPoint3& GetIntersectionPtStart() const
		{
			return m_intersection_start_pt;
		}
		/**
	    * @brief Returns the ending point of the intersection segment.
	    *
	    * @return Constant reference to the ending point.
	    */
		inline const CPoint3& GetIntersectionPtEnd() const
		{
			return m_intersection_end_pt;
		}
		/**
	    * @brief Returns the processing status of the intersection segment.
	    *
	    * Indicates whether the segment has already been used during
	    * the reconstruction of polylines after slicing.
	    *
	    * @return True if the segment has already been used, false otherwise.
	    */
		inline bool GetUsed() const
		{
			return m_segment_is_used;
		}
		/**
		* @brief Sets the processing status of the intersection segment.
		*
		* @param f_segment_is_used True if the segment has been processed,
		*                          false otherwise.
		*/
		inline void SetUsed(const bool f_segment_is_used)
		{
			m_segment_is_used = f_segment_is_used;
		}

		/**
		* @brief Returns the identifier of the mesh edge associated with
		*        the starting point of the segment.
		*
		* @return Identifier of the starting mesh edge.
		*/
		inline uint64_t GetStartEdgeIndex() const
		{
			return m_start_edge_id;
		}
		/**
		* @brief Returns the identifier of the mesh edge associated with
		*        the ending point of the segment.
		*
		* @return Identifier of the ending mesh edge.
		*/
		inline uint64_t GetEndEdgeIndex() const
		{
			return m_end_edge_id;
		}

	private:

		/**
		* @brief Starting point of the oriented intersection segment.
		*
		* Stores the 3D coordinates of the first endpoint of the segment.
		*/
		CPoint3 m_intersection_start_pt;
		/**
		* @brief Ending point of the oriented intersection segment.
		*
		* Stores the 3D coordinates of the second endpoint of the segment.
		*/
		CPoint3 m_intersection_end_pt;
		/**
		* @brief Identifier of the mesh edge associated with the starting point.
		*
		* Identifies the original mesh edge on which the starting intersection
		* point lies.
		*
		* If the starting point corresponds to an existing mesh vertex rather
		* than an intersection occurring inside a mesh edge, this identifier
		* is set to -1.
		*/
		uint64_t m_start_edge_id;
		/**
		* @brief Identifier of the mesh edge associated with the ending point.
		*
		* Identifies the original mesh edge on which the ending intersection
		* point lies.
		*
		* If the ending point corresponds to an existing mesh vertex rather
		* than an intersection occurring inside a mesh edge, this identifier
		* is set to -1.
		*/
		uint64_t m_end_edge_id;
		/**
		* @brief Identifier of the mesh face (triangle) that generated the segment.
		*
		* Stores the index of the original mesh face responsible for generating
		* this intersection segment.
		*/
		int m_face_id;
		/**
		* @brief Indicates whether the intersection segment has already been used.
		*
		* This flag is used during the reconstruction of polylines from the
		* intersection segments produced by mesh slicing.
		*
		* It is initialized to false and set to true once the segment has
		* been incorporated into a reconstructed polyline.
		*/
		bool m_segment_is_used;
	};

	/**
	* @brief Computes the signed distance of a point from a plane.
	*
	* Calculates the signed distance-like value between a point and a plane
	* using the plane normal and a reference point belonging to the plane.
	*
	* The result is positive when the point lies in the direction of the
	* plane normal, negative when it lies in the opposite direction, and
	* zero when the point belongs to the plane.
	*
	* @param f_plane_normal Normal vector of the plane.
	* @param f_any_point_on_plane Any point belonging to the plane.
	* @param f_pt Point whose signed distance from the plane is calculated.
	*
	* @return Signed distance-like value of the point from the plane.
	*
	* @note The plane normal is not normalized by this function.
	*       Consequently, the returned value represents the signed
	*       perpendicular distance only when the normal is a unit vector.
	*/
	inline static double PlanePointDistance(const CPoint3& f_plane_normal,
										    const CPoint3& f_any_point_on_plane,
										    const CPoint3& f_pt)
	{
		const double l_d_coeff = Dot(f_plane_normal, f_any_point_on_plane);
		return Dot(f_plane_normal, f_pt) - l_d_coeff;
	}

	/**
	* @brief Computes the intersection between a plane and a line segment.
	*
	* Determines whether a finite line segment intersects a plane and,
	* if an intersection exists, calculates the corresponding intersection
	* point.
	*
	* The intersection is calculated parametrically using the segment
	* direction and the plane equation:
	*
	*   Dot(N, X) = Dot(N, P)
	*
	* where N is the plane normal, P is a point belonging to the plane,
	* and X is a point on the plane.
	*
	* The function rejects degenerate segments and segments that are
	* parallel or nearly parallel to the plane. It also verifies that
	* the calculated intersection parameter lies within the segment
	* bounds, allowing a small numerical tolerance.
	*
	* @param f_plane_normal Normal vector of the plane.
	* @param f_any_point_on_plane Any point belonging to the plane.
	* @param f_seg_start Starting point of the line segment.
	* @param f_seg_end Ending point of the line segment.
	* @param f_int_pt Output parameter receiving the intersection point
	*                 when an intersection is found.
	 @param f_parallelism_threshold Threshold used to determine whether
	*                                the segment direction is parallel
	*                                or nearly parallel to the plane.
	*                                The default value is 1.0e-10.
	*
	* @return True if the finite segment intersects the plane, false otherwise.
	*
	* @note The intersection point is calculated only when the function
	*       returns true. If the segment is parallel to the plane or
	*       degenerate, the output parameter remains unchanged.
	*
	* @note The plane normal is expected to be non-zero. It does not
	*       need to be normalized because the segment direction is
	*       normalized before calculating the intersection parameter.
	*
	* @note Intersections occurring at the segment endpoints are accepted.
	*
	* @par Reference
	* Philip J. Schneider and David H. Eberly,
	* "Geometric Tools for Computer Graphics",
	* Morgan Kaufmann, 2003, p. 484.
	*
	* The segment-plane intersection calculation is based on the
	* method described in the referenced work. The implementation
	* additionally incorporates numerical tolerances, degenerate
	* segment detection and segment-bound clamping.
	*/
	inline static bool PlaneSegIntersection(const CPoint3& f_plane_normal,
										    const CPoint3& f_any_point_on_plane,
										    const CPoint3& f_seg_start,
										    const CPoint3& f_seg_end,
										    CPoint3& f_int_pt,
										    const double f_parallelism_threshold = 1.0e-10)
	{
		const double l_d = Dot(f_plane_normal, f_any_point_on_plane);
		double l_segm_len = Distance(f_seg_end, f_seg_start);
		if (l_segm_len < 1.0e-12)
		{
			return false;
		}
		const CPoint3 l_segm_dir = (f_seg_end - f_seg_start).GetUnitVector();

		double l_denominator = Dot(l_segm_dir, f_plane_normal);
		if (fabs(l_denominator) < f_parallelism_threshold)
		{
			return false;
		}

		double l_param = (l_d - Dot(f_plane_normal, f_seg_start)) / l_denominator;
		if (l_param < -1.0e-10 || l_param > l_segm_len + 1.0e-10)
		{
			return false;
		}
		l_param = std::clamp(l_param, 0.0, l_segm_len);

		f_int_pt = f_seg_start + l_param * l_segm_dir;

		return true;
	}

	/**
	* @brief Generates a unique identifier for an undirected mesh edge.
	*
	* Creates a 64-bit identifier from the indices of the two vertices
	* defining a mesh edge.
	*
	* The vertex indices are sorted before being encoded, ensuring that
	* the same identifier is generated regardless of the order in which
	* the vertices are provided.
	*
	* The smaller vertex index is stored in the upper 32 bits, while
	* the larger vertex index is stored in the lower 32 bits.
	*
	* This allows an edge to be identified independently of its
	* orientation.
	*
	* @param f_v1 Index of the first vertex defining the edge.
	* @param f_v2 Index of the second vertex defining the edge.
	*
	* @return 64-bit identifier uniquely representing the undirected edge.
	*
	* @note Vertex indices are assumed to fit within 32 bits.
	*/
	inline static std::uint64_t MakeEdgeId(unsigned int f_v1, unsigned int f_v2)
	{
		const unsigned int l_a = std::min(f_v1, f_v2);
		const unsigned int l_b = std::max(f_v1, f_v2);

		return (static_cast<std::uint64_t>(l_a) << 32) | static_cast<std::uint64_t>(l_b);
	}

	/**
	* @brief Appends a 3D point to a flattened profile coordinate vector.
	*
	* Inserts the X, Y and Z coordinates of the supplied point at the
	* end of the profile coordinate vector.
	*
	* Coordinates are stored consecutively as X, Y, Z triplets, allowing
	* multiple 3D points to be represented using a single vector of doubles.
	*
	* @param f_pt 3D point to append to the profile.
	* @param f_prf_coords Vector receiving the point coordinates.
	*
	* @note The input vector is modified by appending exactly three values.
	*/
	inline static void PushPointInProfile(const CPoint3& f_pt, std::vector<double>& f_prf_coords)
	{
		f_prf_coords.emplace_back(f_pt.m_x);
		f_prf_coords.emplace_back(f_pt.m_y);
		f_prf_coords.emplace_back(f_pt.m_z);
	}

	/**
	* @brief Returns the first point of a profile.
	*
	* Extracts the first three coordinates from a flattened profile
	* coordinate array and returns them as a 3D point.
	*
	* @param f_prf_coords Vector containing the profile coordinates,
	*                     stored as consecutive X, Y, Z triplets.
	*
	* @return The first point of the profile. Returns a default-constructed
	*         CPoint3 if the input vector contains fewer than three coordinates.
	*/
	inline static CPoint3 GetProfileFront(std::vector<double> f_prf_coords)
	{
		if (f_prf_coords.size() < 3)
		{
			return CPoint3();
		}
		return CPoint3(f_prf_coords[0], f_prf_coords[1], f_prf_coords[2]);
	}

	/**
	* @brief Returns the last point of a profile
	*
	* Extracts the last three coordinates from a flattened profile
	* coordinate array and returns them as a 3D point
	*
	* @param f_prf_coords Vector containing the profile coordinates,
	*                     stored as consecutive X, Y, Z triplets
	*
	* @return The last point of the profile. Returns a default-constructed
	*         CPoint3 if the input vector contains fewer than three coordinates
	*/
	inline static CPoint3 GetProfileBack(std::vector<double> f_prf_coords)
	{
		if (f_prf_coords.size() < 3)
		{
			return CPoint3();
		}
		const size_t l_size = f_prf_coords.size();
		return CPoint3(f_prf_coords[l_size - 3], f_prf_coords[l_size - 2], f_prf_coords[l_size - 1]);
	}

	/**
	* @brief Determines the orientation and angle of three consecutive points.
	*
	* Computes the orientation of the triplet (P, Q, R) in the XY plane
	* and calculates the angle formed at point Q. The orientation is determined
	* using the signed 2D cross-product of the vectors associated with the three 
	* points.
	*
	* The orientation classification is:
	* - 0: Collinear points.
	* - 1: Clockwise orientation (concave vertex).
	* - 2: Counter-clockwise orientation (convex vertex).
	*
	* The returned angle is expressed in degrees. For counter-clockwise
	* orientations, the angle is converted to its reflex-angle equivalent
	* by subtracting the internal angle from 360 degrees.
	*
	* @param f_p First point of the triplet.
	* @param f_q Middle point of the triplet, where the angle is evaluated.
	* @param f_r Third point of the triplet.
	* @param f_orientation Output parameter containing the orientation
	*                      classification (0, 1 or 2).
	*
	* @return Angle at point Q, expressed in degrees.
	*
	* @note Collinearity is determined using a tolerance of 1.0e-10
	*       applied to the signed orientation value.
	*/
	inline static double TripletOrientationAndAngle(const CPoint3& f_p,
												    const CPoint3& f_q,
												    const CPoint3& f_r,
												    int& f_orientation)
	{
		double l_val = (f_q.m_y - f_p.m_y) * (f_r.m_x - f_q.m_x) - (f_q.m_x - f_p.m_x) * (f_r.m_y - f_q.m_y);

		double l_angle = 0.0;
		l_angle = acos(Dot((f_p - f_q).GetUnitVector(), (f_r - f_q).GetUnitVector()));
		l_angle = l_angle * constants::RAD_TO_DEG;

		if (fabs(l_val) < 1.0e-10)
		{
			f_orientation = 0;
			l_val = 0.0;
			return l_angle;
		}
		if (l_val > 0.0)
		{
			f_orientation = 1;
		}
		else
		{
			f_orientation = 2;
			l_angle = (-l_angle) + 360.0;
		}

		return l_angle;
	}

	/**
	* @brief Finds the point of a bounding box that is farthest from a plane.
	*
	* Determines the bounding-box vertex that maximizes the dot product
	* between the plane normal and the vertex position.
	*
	* The selected vertex is the support point of the axis-aligned bounding
	* box in the direction of the plane normal. Its signed distance from
	* the plane is then calculated.
	*
	* @param f_plane_normal Normal vector of the plane.
	* @param f_any_point_on_plane Any point belonging to the plane.
	* @param f_bbox_ll_vert Lower-left (minimum-coordinate) vertex of the
	*                       axis-aligned bounding box.
	* @param f_box_ur_vert Upper-right (maximum-coordinate) vertex of the
	*                      axis-aligned bounding box.
	* @param f_distance Output parameter containing the signed distance
	*                   between the plane and the selected bounding-box vertex.
	*
	* @return Bounding-box vertex maximizing the dot product with the
	*         plane normal.
	*
	* @note The bounding box is assumed to be axis-aligned, with its
	*       lower and upper vertices defined by the component-wise
	*       minimum and maximum coordinates, respectively.
	*/
	inline static CPoint3 GetBoundingBoxFarthestPoint(const CPoint3& f_plane_normal,
												      const CPoint3& f_any_point_on_plane,
												      const CPoint3& f_bbox_ll_vert,
												      const CPoint3& f_box_ur_vert,
												      double& f_distance)
	{
		CPoint3 l_farthest_point;

		if (f_plane_normal.m_x >= 0.0)
		{
			l_farthest_point.m_x = f_box_ur_vert.m_x;
		}
		else
		{
			l_farthest_point.m_x = f_bbox_ll_vert.m_x;
		}
		if (f_plane_normal.m_y >= 0.0)
		{
			l_farthest_point.m_y = f_box_ur_vert.m_y;
		}
		else
		{
			l_farthest_point.m_y = f_bbox_ll_vert.m_y;
		}
		if (f_plane_normal.m_z >= 0.0)
		{
			l_farthest_point.m_z = f_box_ur_vert.m_z;
		}
		else
		{
			l_farthest_point.m_z = f_bbox_ll_vert.m_z;
		}

		f_distance = PlanePointDistance(f_plane_normal, f_any_point_on_plane, l_farthest_point);

		return l_farthest_point;
	}

	/**
	* @brief Finds the point of a bounding box that is closest to a plane.
	*
	* Determines the bounding-box vertex that minimizes the dot product
	* between the plane normal and the vertex position.
	*
	* The selected vertex is the support point of the axis-aligned bounding
	* box in the direction opposite to the plane normal. Its signed distance
	* from the plane is then calculated.
	*
	* @param f_plane_normal Normal vector of the plane.
	* @param f_any_point_on_plane Any point belonging to the plane.
	* @param f_bbox_ll_vert Lower-left (minimum-coordinate) vertex of the
	*                       axis-aligned bounding box.
	* @param f_box_ur_vert Upper-right (maximum-coordinate) vertex of the
	*                      axis-aligned bounding box.
	* @param f_distance Output parameter containing the signed distance
	*                   between the plane and the selected bounding-box vertex.
	*
	* @return Bounding-box vertex minimizing the dot product with the
	*         plane normal.
	*
	* @note The bounding box is assumed to be axis-aligned, with its
	*       lower and upper vertices defined by the component-wise
	*       minimum and maximum coordinates, respectively.
	*/
	inline static CPoint3 GetBoundingBoxClosestPoint(const CPoint3& f_plane_normal,
												     const CPoint3& f_any_point_on_plane,
												     const CPoint3& f_bbox_ll_vert,
												     const CPoint3& f_box_ur_vert,
												     double& f_distance)
	{
		CPoint3 l_closest_point;

		if (f_plane_normal.m_x >= 0.0)
		{
			l_closest_point.m_x = f_bbox_ll_vert.m_x;
		}
		else
		{
			l_closest_point.m_x = f_box_ur_vert.m_x;
		}

		if (f_plane_normal.m_y >= 0.0)
		{
			l_closest_point.m_y = f_bbox_ll_vert.m_y;
		}
		else
		{
			l_closest_point.m_y = f_box_ur_vert.m_y;
		}

		if (f_plane_normal.m_z >= 0.0)
		{
			l_closest_point.m_z = f_bbox_ll_vert.m_z;
		}
		else
		{
			l_closest_point.m_z = f_box_ur_vert.m_z;
		}

		f_distance = PlanePointDistance(f_plane_normal, f_any_point_on_plane, l_closest_point);

		return l_closest_point;
	}

	/**
	* @brief Computes the axis-aligned bounding box of a mesh.
	*
	* Iterates through the coordinates of all mesh vertices and determines
	* the minimum and maximum values along the X, Y and Z axes.
	*
	* The input coordinates are expected to be stored as consecutive
	* X, Y, Z triplets.
	*
	* @param f_tri_vertexes_coords Flattened vector containing the coordinates
	*                              of the mesh vertices.
	* @param f_ll_vert Output parameter receiving the lower bounding-box
	*                  vertex, containing the minimum coordinate along
	*                  each axis.
	* @param f_ur_vert Output parameter receiving the upper bounding-box
	*                  vertex, containing the maximum coordinate along
	*                  each axis.
	*
	* @note The input vector must contain at least one complete vertex
	*       and its size must be a multiple of three.
	*/
	inline static void GetMeshBoundingBox(const std::vector<double>& f_tri_vertexes_coords,
									      CPoint3& f_ll_vert,
									      CPoint3& f_ur_vert)
	{
		f_ll_vert = CPoint3(f_tri_vertexes_coords[0], f_tri_vertexes_coords[1], f_tri_vertexes_coords[2]);
		f_ur_vert = f_ll_vert;
		int l_coord_id = 0;
		for (unsigned int i = 0; i < f_tri_vertexes_coords.size(); i++)
		{
			l_coord_id += 1;
			if (l_coord_id == 3)
			{
				l_coord_id = 0;
				if (f_tri_vertexes_coords[i - 2] > f_ur_vert.m_x)
				{
					f_ur_vert.m_x = f_tri_vertexes_coords[i - 2];
				}
				if (f_tri_vertexes_coords[i - 1] > f_ur_vert.m_y)
				{
					f_ur_vert.m_y = f_tri_vertexes_coords[i - 1];
				}
				if (f_tri_vertexes_coords[i] > f_ur_vert.m_z)
				{
					f_ur_vert.m_z = f_tri_vertexes_coords[i];
				}
				if (f_tri_vertexes_coords[i - 2] < f_ll_vert.m_x)
				{
					f_ll_vert.m_x = f_tri_vertexes_coords[i - 2];
				}
				if (f_tri_vertexes_coords[i - 1] < f_ll_vert.m_y)
				{
					f_ll_vert.m_y = f_tri_vertexes_coords[i - 1];
				}
				if (f_tri_vertexes_coords[i] < f_ll_vert.m_z)
				{
					f_ll_vert.m_z = f_tri_vertexes_coords[i];
				}
			}
		}
	}

	/**
	* @brief Computes the axis-aligned bounding box of a triangle.
	*
	* Determines the minimum and maximum coordinates of the three
	* triangle vertices independently along the X, Y and Z axes.
	*
	* @param f_vert_1 First vertex of the triangle.
	* @param f_vert_2 Second vertex of the triangle.
	* @param f_vert_3 Third vertex of the triangle.
	* @param f_tri_ll_vert Output parameter receiving the lower bounding-box
	*                      vertex, containing the minimum coordinate
	*                      along each axis.
	* @param f_tri_ur_vert Output parameter receiving the upper bounding-box
	*                      vertex, containing the maximum coordinate
	*                      along each axis.
	*
	* @note The resulting bounding box is axis-aligned and represents
	*       the smallest such box enclosing the three triangle vertices.
	*/
	inline static void GetTriangleBBox(const CPoint3& f_vert_1,
									   const CPoint3& f_vert_2,
									   const CPoint3& f_vert_3,
									   CPoint3& f_tri_ll_vert,
									   CPoint3& f_tri_ur_vert)
	{
		f_tri_ll_vert.m_x = std::min(f_vert_1.m_x, std::min(f_vert_2.m_x, f_vert_3.m_x));
		f_tri_ll_vert.m_y = std::min(f_vert_1.m_y, std::min(f_vert_2.m_y, f_vert_3.m_y));
		f_tri_ll_vert.m_z = std::min(f_vert_1.m_z, std::min(f_vert_2.m_z, f_vert_3.m_z));

		f_tri_ur_vert.m_x = std::max(f_vert_1.m_x, std::max(f_vert_2.m_x, f_vert_3.m_x));
		f_tri_ur_vert.m_y = std::max(f_vert_1.m_y, std::max(f_vert_2.m_y, f_vert_3.m_y));
		f_tri_ur_vert.m_z = std::max(f_vert_1.m_z, std::max(f_vert_2.m_z, f_vert_3.m_z));
	}

	/**
	* @brief Finalizes and stores a reconstructed slicing profile.
	*
	* Checks whether the supplied profile contains enough coordinates
	* to represent at least three 3D points. If valid, appends the first
	* point to the end of the coordinate vector, explicitly closing
	* the profile.
	*
	* The closed profile is then added to the collection of profiles
	* associated with the current slicing plane.
	*
	* If debug mode is enabled, the profile is drawn, an informational
	* message is logged, and the graphics drawing state is reset.
	*
	* @param f_prf_coords Coordinate vector representing the profile,
	*                     stored as consecutive X, Y, Z triplets.
	*                     The vector is modified by appending the first
	*                     point to close the profile.
	* @param f_prfs_current_plane Collection of reconstructed profiles
	*                             associated with the current slicing plane.
	*                             The finalized profile is appended to this
	*                             collection.
	* @param f_graph_ptr Pointer to the graphics interface used for
	*                    visualization and logging when debug mode is enabled.
	* @param f_debug Enables or disables debug visualization and logging.
	*
	* @note Profiles containing fewer than nine coordinates are ignored.
	*       The input profile is modified in place.
	*/
	inline static void DispatchSlicingProfile(std::vector<double>& f_prf_coords,
										      std::vector<std::vector<double>>& f_prfs_current_plane
											  )
	{
		if (f_prf_coords.size() < 9)
		{
			return;
		}
		f_prf_coords.emplace_back(f_prf_coords[0]);
		f_prf_coords.emplace_back(f_prf_coords[1]);
		f_prf_coords.emplace_back(f_prf_coords[2]);

		f_prfs_current_plane.emplace_back(f_prf_coords);
	}

	/**
	* @brief Computes the intersection between a triangle and a plane.
	*
	* Determines whether a triangle intersects a plane and, if an
	* intersection exists, calculates the two endpoints of the resulting
	* intersection segment.
	*
	* The position of each triangle vertex relative to the plane is
	* determined by evaluating the signed dot product between the plane
	* normal and the vector connecting a point on the plane to the vertex.
	*
	* A positive value indicates that the vertex lies in the direction
	* of the plane normal, a negative value indicates that it lies in
	* the opposite direction, and a value within the specified tolerance
	* is considered to lie on the plane.
	*
	* The function handles the different possible configurations of
	* triangle-plane intersection:
	*
	* - Three vertices on the same side of the plane: no intersection.
	* - Three vertices on the plane: coplanar triangle, ignored.
	* - Two vertices on one side and the third on the opposite side:
	*   the intersection segment is computed between the corresponding
	*   triangle edges.
	* - One vertex on the plane and the other two on the same side:
	*   the intersection is a single point and is discarded.
	* - One vertex on the plane and the other two on opposite sides:
	*   the intersection segment connects the on-plane vertex to the
	*   intersection point on the opposite edge.
	* - Two vertices on the plane: the corresponding triangle edge lies
	*   on the plane. The edge is either accepted or discarded according
	*   to the position of the third vertex.
	*
	* Special handling of vertices and edges lying on the plane is
	* necessary to avoid generating duplicate intersection segments
	* when processing adjacent triangles of a manifold mesh.
	*
	* In particular, when an entire triangle edge lies on the plane,
	* the intersection is considered only if the third vertex lies
	* on the negative side of the plane. This convention is consistent
	* with the slicing strategy adopted by this implementation.
	*
	* The resulting segment endpoints are ordered according to the
	* orientation of the originating triangle. Assuming consistently
	* oriented triangle vertices (counter-clockwise when viewed from
	* outside the mesh), the triangle normal points to the right-hand
	* side of the oriented intersection segment.
	*
	* The edge identifiers associated with the intersection endpoints
	* are also returned. When an endpoint corresponds to an existing
	* triangle vertex rather than an intersection occurring inside a
	* triangle edge, its identifier is set to UINT64_MAX.
	*
	* @param f_tri_vert_1 First vertex of the triangle.
	* @param f_tri_vert_2 Second vertex of the triangle.
	* @param f_tri_vert_3 Third vertex of the triangle.
	* @param f_edge_1_id Identifier of the triangle edge connecting
	*                    vertices 1 and 2.
	* @param f_edge_2_id Identifier of the triangle edge connecting
	*                    vertices 2 and 3.
	* @param f_edge_3_id Identifier of the triangle edge connecting
	*                    vertices 3 and 1.
	* @param f_plane_normal Normal vector of the slicing plane.
	* @param f_any_point_on_plane Any point belonging to the slicing plane.
	* @param f_int_pt_on_edge_start Output parameter receiving the starting
	*                               point of the intersection segment.
	* @param f_int_pt_on_edge_end Output parameter receiving the ending
	*                             point of the intersection segment.
	* @param f_edge_index_start Output parameter receiving the identifier
	*                           of the triangle edge associated with the
	*                           starting intersection point. Set to
	*                           UINT64_MAX if the point corresponds to
	*                           an existing triangle vertex.
	* @param f_edge_index_end Output parameter receiving the identifier
	*                         of the triangle edge associated with the
	*                         ending intersection point. Set to
	*                         UINT64_MAX if the point corresponds to
	*                         an existing triangle vertex.
	* @param f_parallelism_threshold Tolerance used to classify triangle
	*                                vertices as lying on the slicing plane.
	*                                Values whose absolute signed dot product
	*                                does not exceed this threshold are
	*                                considered zero.
	*
	* @return True if a valid intersection segment is generated,
	*         false otherwise.
	*
	* @note The function modifies the output parameters only when a valid
	*       intersection configuration is identified. However, the
	*       intersection calculations delegated to PlaneSegIntersection()
	*       are not checked for their return values.
	*
	* @note A triangle entirely contained in the slicing plane is ignored.
	*
	* @note An intersection consisting of a single point is discarded,
	*       since zero-length segments are not useful for profile
	*       reconstruction.
	*
	* @note When an entire triangle edge lies on the slicing plane,
	*       the intersection is accepted only if the third vertex lies
	*       on the negative side of the plane. Reversing this convention
	*       would change which boundary layer is included when a slicing
	*       plane coincides with the bottom or top of the model.
	*
	* @note The plane normal does not need to be normalized, provided
	*       the tolerance is interpreted consistently with its magnitude.
	*
	* @par Reference
	* Philip J. Schneider and David H. Eberly,
	* "Geometric Tools for Computer Graphics",
	* Morgan Kaufmann, 2003, Chapter 11, Section 11.5,
	* Figures 11.25 and 11.26, pp. 534-535.
	*
	* The general triangle-plane intersection configurations are based
	* on the geometric treatment presented in the referenced work.
	* The handling of coplanar edges, vertex intersections and boundary
	* inclusion follows the specific conventions adopted by this
	* mesh-slicing implementation.
	*/
	inline static bool PlaneTriangleIntersect(const CPoint3& f_tri_vert_1,
										      const CPoint3& f_tri_vert_2,
										      const CPoint3& f_tri_vert_3,
										      const uint64_t f_edge_1_id,
										      const uint64_t f_edge_2_id,
										      const uint64_t f_edge_3_id,
										      const CPoint3& f_plane_normal,
										      const CPoint3& f_any_point_on_plane,
										      CPoint3& f_int_pt_on_edge_start,
										      CPoint3& f_int_pt_on_edge_end,
										      uint64_t& f_edge_index_start,
										      uint64_t& f_edge_index_end,
										      const double f_parallelism_threshold)
	{


		double l_dot_1 = 0.0;
		double l_dot_2 = 0.0;
		double l_dot_3 = 0.0;
		l_dot_1 = Dot(f_plane_normal, (f_tri_vert_1 - f_any_point_on_plane));
		l_dot_2 = Dot(f_plane_normal, (f_tri_vert_2 - f_any_point_on_plane));
		l_dot_3 = Dot(f_plane_normal, (f_tri_vert_3 - f_any_point_on_plane));

		if (fabs(l_dot_1) <= f_parallelism_threshold)
		{
			l_dot_1 = 0.0;
		}
		if (fabs(l_dot_2) <= f_parallelism_threshold)
		{
			l_dot_2 = 0.0;
		}
		if (fabs(l_dot_3) <= f_parallelism_threshold)
		{
			l_dot_3 = 0.0;
		}

		if ((l_dot_1 > 0) && (l_dot_2 > 0) && (l_dot_3 > 0))
		{
			return false;
		}
		else if ((l_dot_1 < 0) && (l_dot_2 < 0) && (l_dot_3 < 0))
		{
			return false;
		}
		if ((fabs(l_dot_1) + fabs(l_dot_2) + fabs(l_dot_3)) < 1.0e-10)
		{
			return false;
		}

		if ((l_dot_1 > 0) && (l_dot_2 > 0) && (l_dot_3 < 0))
		{
			PlaneSegIntersection(f_plane_normal, f_any_point_on_plane, f_tri_vert_2, f_tri_vert_3, f_int_pt_on_edge_start);
			PlaneSegIntersection(f_plane_normal, f_any_point_on_plane, f_tri_vert_1, f_tri_vert_3, f_int_pt_on_edge_end);
			f_edge_index_start = f_edge_2_id;
			f_edge_index_end = f_edge_3_id;

			return true;
		}
		if ((l_dot_1 < 0) && (l_dot_2 < 0) && (l_dot_3 > 0))
		{
			PlaneSegIntersection(f_plane_normal, f_any_point_on_plane, f_tri_vert_1, f_tri_vert_3, f_int_pt_on_edge_start);
			PlaneSegIntersection(f_plane_normal, f_any_point_on_plane, f_tri_vert_2, f_tri_vert_3, f_int_pt_on_edge_end);
			f_edge_index_start = f_edge_3_id;
			f_edge_index_end = f_edge_2_id;

			return true;
		}
		if ((l_dot_2 > 0) && (l_dot_3 > 0) && (l_dot_1 < 0))
		{
			PlaneSegIntersection(f_plane_normal, f_any_point_on_plane, f_tri_vert_3, f_tri_vert_1, f_int_pt_on_edge_start);
			PlaneSegIntersection(f_plane_normal, f_any_point_on_plane, f_tri_vert_2, f_tri_vert_1, f_int_pt_on_edge_end);
			f_edge_index_start = f_edge_3_id;
			f_edge_index_end = f_edge_1_id;

			return true;
		}
		if ((l_dot_2 < 0) && (l_dot_3 < 0) && (l_dot_1 > 0))
		{
			PlaneSegIntersection(f_plane_normal, f_any_point_on_plane, f_tri_vert_2, f_tri_vert_1, f_int_pt_on_edge_start);
			PlaneSegIntersection(f_plane_normal, f_any_point_on_plane, f_tri_vert_3, f_tri_vert_1, f_int_pt_on_edge_end);
			f_edge_index_start = f_edge_1_id;
			f_edge_index_end = f_edge_3_id;

			return true;
		}
		if ((l_dot_1 > 0) && (l_dot_3 > 0) && (l_dot_2 < 0))
		{
			PlaneSegIntersection(f_plane_normal, f_any_point_on_plane, f_tri_vert_1, f_tri_vert_2, f_int_pt_on_edge_start);
			PlaneSegIntersection(f_plane_normal, f_any_point_on_plane, f_tri_vert_3, f_tri_vert_2, f_int_pt_on_edge_end);
			f_edge_index_start = f_edge_1_id;
			f_edge_index_end = f_edge_2_id;

			return true;
		}
		if ((l_dot_1 < 0) && (l_dot_3 < 0) && (l_dot_2 > 0))
		{
			PlaneSegIntersection(f_plane_normal, f_any_point_on_plane, f_tri_vert_3, f_tri_vert_2, f_int_pt_on_edge_start);
			PlaneSegIntersection(f_plane_normal, f_any_point_on_plane, f_tri_vert_1, f_tri_vert_2, f_int_pt_on_edge_end);
			f_edge_index_start = f_edge_2_id;
			f_edge_index_end = f_edge_1_id;

			return true;
		}

		if ((l_dot_1 == 0) && (((l_dot_2 > 0) && (l_dot_3 > 0)) || ((l_dot_2 < 0) && (l_dot_3 < 0))))
		{

			return false;
		}
		if ((l_dot_2 == 0) && (((l_dot_1 > 0) && (l_dot_3 > 0)) || ((l_dot_1 < 0) && (l_dot_3 < 0))))
		{

			return false;
		}
		if ((l_dot_3 == 0) && (((l_dot_2 > 0) && (l_dot_1 > 0)) || ((l_dot_2 < 0) && (l_dot_1 < 0))))
		{

			return false;
		}

		if ((l_dot_1 == 0) && ((l_dot_2 > 0) && (l_dot_3 < 0)))
		{
			PlaneSegIntersection(f_plane_normal, f_any_point_on_plane, f_tri_vert_3, f_tri_vert_2, f_int_pt_on_edge_start);
			f_int_pt_on_edge_end = f_tri_vert_1;

			f_edge_index_start = f_edge_2_id;
			f_edge_index_end = UINT64_MAX;

			return true;
		}
		if ((l_dot_1 == 0) && ((l_dot_2 < 0) && (l_dot_3 > 0)))
		{
			f_int_pt_on_edge_start = f_tri_vert_1;
			PlaneSegIntersection(f_plane_normal, f_any_point_on_plane, f_tri_vert_3, f_tri_vert_2, f_int_pt_on_edge_end);

			f_edge_index_start = UINT64_MAX;
			f_edge_index_end = f_edge_2_id;

			return true;
		}
		if ((l_dot_2 == 0) && ((l_dot_1 > 0) && (l_dot_3 < 0)))
		{
			f_int_pt_on_edge_start = f_tri_vert_2;
			PlaneSegIntersection(f_plane_normal, f_any_point_on_plane, f_tri_vert_1, f_tri_vert_3, f_int_pt_on_edge_end);

			f_edge_index_start = UINT64_MAX;
			f_edge_index_end = f_edge_3_id;

			return true;
		}
		if ((l_dot_2 == 0) && ((l_dot_1 < 0) && (l_dot_3 > 0)))
		{
			PlaneSegIntersection(f_plane_normal, f_any_point_on_plane, f_tri_vert_1, f_tri_vert_3, f_int_pt_on_edge_start);
			f_int_pt_on_edge_end = f_tri_vert_2;

			f_edge_index_start = f_edge_3_id;
			f_edge_index_end = UINT64_MAX;

			return true;
		}
		if ((l_dot_3 == 0) && ((l_dot_1 > 0) && (l_dot_2 < 0)))
		{
			PlaneSegIntersection(f_plane_normal, f_any_point_on_plane, f_tri_vert_1, f_tri_vert_2, f_int_pt_on_edge_start);
			f_int_pt_on_edge_end = f_tri_vert_3;

			f_edge_index_start = f_edge_1_id;
			f_edge_index_end = UINT64_MAX;

			return true;
		}
		if ((l_dot_3 == 0) && ((l_dot_1 < 0) && (l_dot_2 > 0)))
		{
			f_int_pt_on_edge_start = f_tri_vert_3;
			PlaneSegIntersection(f_plane_normal, f_any_point_on_plane, f_tri_vert_1, f_tri_vert_2, f_int_pt_on_edge_end);

			f_edge_index_start = UINT64_MAX;
			f_edge_index_end = f_edge_1_id;

			return true;
		}

		if ((l_dot_1 == 0) && (l_dot_2 == 0))
		{
			if (l_dot_3 > 0)
			{
				return false;
			}
			f_int_pt_on_edge_start = f_tri_vert_2;
			f_int_pt_on_edge_end = f_tri_vert_1;

			f_edge_index_start = UINT64_MAX;
			f_edge_index_end = UINT64_MAX;

			return true;
		}
		if ((l_dot_2 == 0) && (l_dot_3 == 0))
		{
			if (l_dot_1 > 0)
			{
				return false;
			}
			f_int_pt_on_edge_start = f_tri_vert_3;
			f_int_pt_on_edge_end = f_tri_vert_2;

			f_edge_index_start = UINT64_MAX;
			f_edge_index_end = UINT64_MAX;

			return true;
		}
		if ((l_dot_1 == 0) && (l_dot_3 == 0))
		{
			if (l_dot_2 > 0)
			{
				return false;
			}
			f_int_pt_on_edge_start = f_tri_vert_1;
			f_int_pt_on_edge_end = f_tri_vert_3;

			f_edge_index_start = UINT64_MAX;
			f_edge_index_end = UINT64_MAX;

			return true;
		}

		return false;

	}

	/**
	* @brief Slices a triangular mesh into planar profiles at regular intervals
	*
	* Performs the slicing of a triangular mesh using a sequence of parallel
	* planes, generating the intersection profiles associated with each
	* slicing plane. The input mesh is represented by a list of vertex coordina-
	* tes and a list of triangle vertex indices. Each slicing plane is defined by
	* a reference point and a normal vector. Successive slicing planes are
	* displaced along the plane normal according to the specified layer
	* height.
	*
	* For each slicing plane, the function:
	* - Identifies the triangles intersected by the plane.
	* - Computes the intersection segments between the triangles and
	*   the slicing plane.
	* - Reconstructs closed profiles by connecting the intersection
	*   segments.
	* - Stores the resulting profiles in the output container.
	*
	* The reconstructed profiles preserve the orientation induced by
	* the original mesh triangle orientation, allowing external contours
	* and internal contours (holes) to be distinguished.
	*
	* @param f_tri_vertexes_coords		Flattened vector containing the coordinates
	*									of the mesh vertices, stored as consecutive
	*									X, Y, Z triplets.
	* @param f_tri_indexes				Vector containing the vertex indices defining
	*									the mesh triangles. Each consecutive group of
	*									three indices identifies one triangle.
	* @param f_layer_height				Distance between two consecutive slicing planes,
	*									measured along the slicing direction.
	* @param f_start_point_x			X coordinate of a reference point belonging
	*									to the first slicing plane.
	* @param f_start_point_y			Y coordinate of a reference point belonging
	*									to the first slicing plane.
	* @param f_start_point_z			Z coordinate of a reference point belonging
	*									to the first slicing plane.
	* @param f_slicing_plane_normal_x	X component of the slicing plane normal.
	* @param f_slicing_plane_normal_y	Y component of the slicing plane normal.
	* @param f_slicing_plane_normal_z	Z component of the slicing plane normal.
	* @param f_resulting_profiles_coords Output container receiving the
	*                                    reconstructed profiles for each
	*                                    slicing plane.
	*                                    The container has three levels:
	*                                    - Outer vector: slicing planes.
	*                                    - Middle vector: profiles generated
	*                                      by the corresponding slicing plane.
	*                                    - Inner vector: flattened coordinates
	*                                      of the profile vertices, stored as
	*                                      consecutive X, Y, Z triplets.
	*
	*                                    Each reconstructed profile is
	*                                    represented by a sequence of 3D
	*                                    coordinates, with the first point
	*                                    repeated at the end to explicitly
	*                                    close the profile.
	*
	* @return	Integer status code indicating the outcome of the slicing
	*			operation.
	*			0 = slicing ok
	*			1 = invalid input parameters
	*			2 = error during slicing
	* @note The slicing planes are parallel and their positions are
	*       determined by the reference point, the slicing direction
	*       and the layer height.
	* @note The orientation of the slicing plane normal determines the
	*       direction in which successive slicing planes are generated.
	* @note The input mesh is expected to be represented by valid triangle
	*       indices referencing existing vertices.
	* @note The interpretation of the resulting profiles as external
	*       contours or internal contours depends on the consistent
	*       orientation of the input mesh triangles.
	* @note Debug visualization requires a valid graphics interface pointer
	*       whenever the selected debugging level activates drawing or
	*       logging operations.
	*/
	inline int SliceTriangularMesh(const std::vector<double>& f_tri_vertexes_coords,
								   const std::vector<unsigned int>& f_tri_indexes,
								   const double f_layer_height,
								   const double f_start_point_x,
								   const double f_start_point_y,
								   const double f_start_point_z,
								   const double f_slicing_plane_normal_x,
								   const double f_slicing_plane_normal_y,
								   const double f_slicing_plane_normal_z,
								   std::vector<std::vector<std::vector<double>>>& f_resulting_profiles_coords
								   )
	{

		if ((f_tri_vertexes_coords.size() < 3) ||
			(f_tri_indexes.size() < 3) ||
			(f_layer_height < 1.0e-12))
		{
			return 1;
		}

		f_resulting_profiles_coords.clear();
		std::vector<std::unordered_map<uint64_t, std::vector<CIntersectionSegment>>> l_segment_intersection_map_planes;

		CPoint3 l_slicing_plane_normal(f_slicing_plane_normal_x, f_slicing_plane_normal_y, f_slicing_plane_normal_z);
		l_slicing_plane_normal.UnitizeInPlace();
		const CPoint3 l_point_on_first_slicing_plane(f_start_point_x, f_start_point_y, f_start_point_z);

		CPoint3 l_ll_vert;
		CPoint3 l_ur_vert;
		GetMeshBoundingBox(f_tri_vertexes_coords, l_ll_vert, l_ur_vert);

		double l_farthest_pt_dist = 0.0;
		CPoint3 l_start_plane_farthest_pt = GetBoundingBoxFarthestPoint(l_slicing_plane_normal,
																		l_point_on_first_slicing_plane,
																		l_ll_vert,
																		l_ur_vert,
																		l_farthest_pt_dist);
		if (l_farthest_pt_dist <= 0.0)
		{
			return 0;
		}

		const double l_slicing_span = PlanePointDistance(l_slicing_plane_normal,
														 l_point_on_first_slicing_plane,
														 l_start_plane_farthest_pt);
		int l_num_division = (int)(l_slicing_span / f_layer_height);
		l_segment_intersection_map_planes.reserve(l_num_division + 2);
		for (int i = 0; i < l_num_division + 2; i++)
		{
			std::unordered_map <uint64_t, std::vector<CIntersectionSegment>> l_intersezioni;
			l_segment_intersection_map_planes.emplace_back(l_intersezioni);
		}

		double l_z_intersezione = 0.0;
		int l_tot_triangles = int(f_tri_indexes.size() / 3);
		for (int l_itriangle = 0; l_itriangle < l_tot_triangles; l_itriangle++)
		{
			unsigned int l_i0 = f_tri_indexes[l_itriangle * 3 + 0];
			unsigned int l_i1 = f_tri_indexes[l_itriangle * 3 + 1];
			unsigned int l_i2 = f_tri_indexes[l_itriangle * 3 + 2];

			const double l_x0 = f_tri_vertexes_coords[l_i0 * 3 + 0];
			const double l_y0 = f_tri_vertexes_coords[l_i0 * 3 + 1];
			const double l_z0 = f_tri_vertexes_coords[l_i0 * 3 + 2];

			const double l_x1 = f_tri_vertexes_coords[l_i1 * 3 + 0];
			const double l_y1 = f_tri_vertexes_coords[l_i1 * 3 + 1];
			const double l_z1 = f_tri_vertexes_coords[l_i1 * 3 + 2];

			const double l_x2 = f_tri_vertexes_coords[l_i2 * 3 + 0];
			const double l_y2 = f_tri_vertexes_coords[l_i2 * 3 + 1];
			const double l_z2 = f_tri_vertexes_coords[l_i2 * 3 + 2];

			const CPoint3 l_vert1(l_x0, l_y0, l_z0);
			const CPoint3 l_vert2(l_x1, l_y1, l_z1);
			const CPoint3 l_vert3(l_x2, l_y2, l_z2);

			const std::uint64_t l_edge1_id = MakeEdgeId(l_i0, l_i1);
			const std::uint64_t l_edge2_id = MakeEdgeId(l_i1, l_i2);
			const std::uint64_t l_edge3_id = MakeEdgeId(l_i2, l_i0);

			CPoint3 l_tri_ll_vert;
			CPoint3 l_tri_ur_vert;
			GetTriangleBBox(l_vert1, l_vert2, l_vert3, l_tri_ll_vert, l_tri_ur_vert);

			double l_closest_dist = 0.0;
			double l_farthest_dist = 0.0;
			CPoint3 l_tri_farthest_pt = GetBoundingBoxFarthestPoint(l_slicing_plane_normal,
																	l_point_on_first_slicing_plane,
																	l_tri_ll_vert,
																	l_tri_ur_vert,
																	l_farthest_dist);
			CPoint3 l_tri_closest_pt = GetBoundingBoxClosestPoint(l_slicing_plane_normal,
																  l_point_on_first_slicing_plane,
																  l_tri_ll_vert,
																  l_tri_ur_vert,
																  l_closest_dist);
			if ((l_closest_dist < 0.0) && (l_farthest_dist < 0.0))
			{
				continue;
			}
			if (l_closest_dist)
			{
				l_closest_dist = 0.0;
			}
			int iplane_start = (int)(l_closest_dist / f_layer_height);
			int iplane_end = (int)(l_farthest_dist / f_layer_height);
			for (int l_kplane = iplane_start; l_kplane < iplane_end + 2; l_kplane++)
			{
				l_z_intersezione = l_kplane * f_layer_height;
				CPoint3 l_point_on_plane = l_point_on_first_slicing_plane + l_slicing_plane_normal * l_z_intersezione;
				CPoint3 l_edge_start_pt;
				CPoint3 l_edge_end_pt;
				uint64_t l_edge_index_start = UINT64_MAX;
				uint64_t l_edge_index_end = UINT64_MAX;
				bool l_ret_int = PlaneTriangleIntersect(l_vert1,
														l_vert2,
														l_vert3,
														l_edge1_id,
														l_edge2_id,
														l_edge3_id,
														l_slicing_plane_normal,
														l_point_on_plane,
														l_edge_start_pt,
														l_edge_end_pt,
														l_edge_index_start,
														l_edge_index_end,
														1.0e-10);
				if (l_ret_int)
				{
					l_segment_intersection_map_planes[l_kplane][l_edge_index_start].push_back(CIntersectionSegment(l_edge_start_pt,
																												   l_edge_end_pt,
																												   l_edge_index_start,
																												   l_edge_index_end,
																												   l_itriangle));
				}
			}
		}
		

		int l_tot_planes = (int)l_segment_intersection_map_planes.size();

		for (int l_iplane = 0; l_iplane < l_tot_planes; l_iplane++)
		{

			std::vector<std::vector<double>> l_prf_list_current_plane;
			if (l_segment_intersection_map_planes[l_iplane].size() == 0)
			{
				continue;
			}

			while (1)
			{
				bool l_seg_found = false;
				uint64_t l_current_ending_index = UINT64_MAX;
				uint64_t l_current_starting_index = UINT64_MAX;
				CIntersectionSegment l_curr_int_seg;

				for (auto l_map_it = l_segment_intersection_map_planes[l_iplane].begin(); l_map_it != l_segment_intersection_map_planes[l_iplane].end(); l_map_it++)
				{
					for (auto l_seg_it = l_map_it->second.begin(); l_seg_it != l_map_it->second.end(); l_seg_it++)
					{
						if (!l_seg_it->GetUsed())
						{
							l_current_ending_index = l_seg_it->GetEndEdgeIndex();
							l_current_starting_index = l_seg_it->GetStartEdgeIndex();
							l_curr_int_seg = *l_seg_it;
							l_seg_it->SetUsed(true);


							l_seg_found = true;
							break;
						}
					}
					if (l_seg_found)
					{
						break;
					}
				}
				if (!l_seg_found)
				{
					break;
				}

				std::vector<double> l_new_profile;
				PushPointInProfile(l_curr_int_seg.GetIntersectionPtStart(), l_new_profile);
				PushPointInProfile(l_curr_int_seg.GetIntersectionPtEnd(), l_new_profile);

				while (1)
				{
					auto l_map_it = l_segment_intersection_map_planes[l_iplane].find(l_current_ending_index);
					if (l_map_it != l_segment_intersection_map_planes[l_iplane].end())
					{
						if (l_current_ending_index == UINT64_MAX)
						{
							std::vector<int> l_candidate_segments_indexes;
							int l_curr_index = 0;
							for (auto l_seg_it = l_map_it->second.begin(); l_seg_it != l_map_it->second.end(); l_seg_it++)
							{
								if ((Distance(l_curr_int_seg.GetIntersectionPtEnd(), l_seg_it->GetIntersectionPtStart()) < 1.0e-10) &&
									(!l_seg_it->GetUsed()))
								{
									l_candidate_segments_indexes.push_back(l_curr_index);
								}

								l_curr_index += 1;
							}
							if (l_candidate_segments_indexes.size() > 1)
							{
								double l_turnleft_metric = -HUGE_VAL;
								int l_selected_seg_id = -1;
								CPoint3 l_p = l_curr_int_seg.GetIntersectionPtStart();
								CPoint3 l_q = l_curr_int_seg.GetIntersectionPtEnd();
								for (auto l_seg_index_it = l_candidate_segments_indexes.begin(); l_seg_index_it != l_candidate_segments_indexes.end(); l_seg_index_it++)
								{
									CPoint3 l_r = l_map_it->second.at(*l_seg_index_it).GetIntersectionPtEnd();
									int l_orientation = 0;
									double l_angle_deg = TripletOrientationAndAngle(l_p, l_q, l_r, l_orientation);
									if (l_angle_deg > l_turnleft_metric)
									{
										l_selected_seg_id = *l_seg_index_it;
										l_turnleft_metric = l_angle_deg;
									}
								}
								if (l_selected_seg_id == -1)
								{
									if (Distance(l_curr_int_seg.GetIntersectionPtEnd(), GetProfileFront(l_new_profile)) < 1.0e-10)
									{
										DispatchSlicingProfile(l_new_profile, 
															   l_prf_list_current_plane 
															   );

										break;
									}
									else
									{
										return 2;
									}
								}
								else
								{

									if (Distance(GetProfileFront(l_new_profile), (l_map_it->second[l_selected_seg_id].GetIntersectionPtEnd())) < 1.0e-10)
									{
										DispatchSlicingProfile(l_new_profile, 
											l_prf_list_current_plane
											);
										l_map_it->second[l_selected_seg_id].SetUsed(true);

										break;
									}
									else
									{
										l_curr_int_seg = l_map_it->second[l_selected_seg_id];
										l_map_it->second[l_selected_seg_id].SetUsed(true);
										l_current_starting_index = l_curr_int_seg.GetStartEdgeIndex();
										l_current_ending_index = l_curr_int_seg.GetEndEdgeIndex();
										PushPointInProfile(l_curr_int_seg.GetIntersectionPtEnd(), l_new_profile);
									}
								}
							}
							else if (l_candidate_segments_indexes.size() == 0)
							{
								if (Distance(l_curr_int_seg.GetIntersectionPtEnd(), GetProfileFront(l_new_profile)) < 1.0e-10)
								{
									DispatchSlicingProfile(l_new_profile, 
														   l_prf_list_current_plane
														   );
									break;
								}
								else
								{
									return 2;
								}
							}
							else if (l_candidate_segments_indexes.size() == 1)
							{
								if (Distance(GetProfileFront(l_new_profile), (l_map_it->second.at(l_candidate_segments_indexes.front()).GetIntersectionPtEnd())) < 1.0e-10)
								{
									DispatchSlicingProfile(l_new_profile, 
														   l_prf_list_current_plane
														   );
									l_map_it->second.at(l_candidate_segments_indexes.front()).SetUsed(true);

									break;
								}
								else
								{
									l_curr_int_seg = l_map_it->second.at(l_candidate_segments_indexes.front());
									l_map_it->second.at(l_candidate_segments_indexes.front()).SetUsed(true);
									l_current_starting_index = l_curr_int_seg.GetStartEdgeIndex();
									l_current_ending_index = l_curr_int_seg.GetEndEdgeIndex();
									PushPointInProfile(l_curr_int_seg.GetIntersectionPtEnd(), l_new_profile);
								}
							}
							else
							{
								
								return 2;
							}
						}
						else if (l_map_it->second.size() > 1)
						{

							double l_turnleft_metric = -HUGE_VAL;
							int l_selected_seg_id = -1;
							int l_seg_id = 0;
							CPoint3 l_p = l_curr_int_seg.GetIntersectionPtStart();
							CPoint3 l_q = l_curr_int_seg.GetIntersectionPtEnd();
							for (auto l_seg_it = l_map_it->second.begin(); l_seg_it != l_map_it->second.end(); l_seg_it++)
							{
								int l_orientation = 0;
								CPoint3 l_r = l_seg_it->GetIntersectionPtEnd();
								double l_angle_deg = TripletOrientationAndAngle(l_p, l_q, l_r, l_orientation);
								if ((l_angle_deg > l_turnleft_metric) && (!l_seg_it->GetUsed()))
								{
									l_selected_seg_id = l_seg_id;
									l_turnleft_metric = l_angle_deg;
								}

								l_seg_id += 1;
							}

							if (l_selected_seg_id == -1)
							{
								if (Distance(l_curr_int_seg.GetIntersectionPtEnd(), GetProfileFront(l_new_profile)) < 1.0e-10)
								{
									DispatchSlicingProfile(l_new_profile, 
														   l_prf_list_current_plane
														   );

									break;
								}
								else
								{
									return 2;
								}
							}
							else
							{

								if (Distance(GetProfileFront(l_new_profile), (l_map_it->second[l_selected_seg_id].GetIntersectionPtEnd())) < 1.0e-10)
								{
									DispatchSlicingProfile(l_new_profile, 
														   l_prf_list_current_plane
														   );
									l_map_it->second[l_selected_seg_id].SetUsed(true);

									break;
								}
								else
								{
									l_curr_int_seg = l_map_it->second[l_selected_seg_id];
									l_map_it->second[l_selected_seg_id].SetUsed(true);
									l_current_starting_index = l_curr_int_seg.GetStartEdgeIndex();
									l_current_ending_index = l_curr_int_seg.GetEndEdgeIndex();
									PushPointInProfile(l_curr_int_seg.GetIntersectionPtEnd(), l_new_profile);
								}
							}
						}
						else
						{

							if (Distance(GetProfileFront(l_new_profile), l_map_it->second.front().GetIntersectionPtEnd()) < 1.0e-10)
							{
								DispatchSlicingProfile(l_new_profile, 
													   l_prf_list_current_plane
													   );
								l_map_it->second.front().SetUsed(true);

								break;
							}
							else
							{
								l_curr_int_seg = l_map_it->second.front();
								l_map_it->second.front().SetUsed(true);
								l_current_starting_index = l_curr_int_seg.GetStartEdgeIndex();
								l_current_ending_index = l_curr_int_seg.GetEndEdgeIndex();
								PushPointInProfile(l_curr_int_seg.GetIntersectionPtEnd(), l_new_profile);
							}
						}
					}
					else
					{
						return 2;
					}
				}
			}

			f_resulting_profiles_coords.emplace_back(l_prf_list_current_plane);
		}

		return 0;
	}
}