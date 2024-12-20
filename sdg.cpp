// standard includes
#include <iostream>
#include <fstream>
#include <cassert>
#include <string>
#include <chrono>
#include <ctime>

// define the kernel
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

//typedef CGAL::Exact_predicates_inexact_constructions_kernel EPIC_Kernel;

// choose the kernel
#include <CGAL/Simple_cartesian.h>
//typedef CGAL::Simple_cartesian<double> SC_Kernel;

#include <CGAL/Lazy_exact_nt.h>
//typedef CGAL::internal::Exact_type_selector<int>::Type Exact_nt;
//typedef CGAL::Simple_cartesian< CGAL::Lazy_exact_nt<Exact_nt> > LZNT_Kernel;
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>

//typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt EK;

//typedef LZNT_Kernel Kernel;

// typedefs for the traits and the algorithm
#include <CGAL/Segment_Delaunay_graph_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_2.h>
#include <CGAL/Segment_Delaunay_graph_filtered_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_storage_traits_with_info_2.h>

//typedef CGAL::Segment_Delaunay_graph_traits_2<Kernel>  Gt;
//typedef CGAL::Segment_Delaunay_graph_2<Gt>             SDG;

//#define SDG_DRAW_DEBUG

typedef CGAL::Simple_cartesian<double>                                            CK;
typedef CK::Point_2                                                               Point_2;
typedef CK::Segment_2                                                             Segment_2;
typedef CGAL::Field_with_sqrt_tag                                                 CM;
typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt               EK;
typedef CGAL::Field_with_sqrt_tag                                                 EM;
typedef CGAL::Simple_cartesian<CGAL::Interval_nt<false> >                         FK;
typedef CGAL::Field_with_sqrt_tag                                                 FM;
typedef CGAL::Segment_Delaunay_graph_filtered_traits_2<CK, CM, EK, EM, FK, FM>    Gt;
//typedef CGAL::Segment_Delaunay_graph_2<Gt>                                        SDG;

typedef CK Kernel;

typedef int Red_blue;

// functor that defines how to convert color info when:
// 1. constructing the storage site of an endpoint of a segment
// 2. a segment site is split into two sub-segments
struct Red_blue_convert_info
{
  typedef Red_blue      Info;
  typedef const Info&   result_type;

  inline
  const Info& operator()(const Info& info0, bool) const {
    // just return the info of the supporting segment
    return info0;
  }

  inline
  const Info& operator()(const Info& info0, const Info& , bool) const {
    // just return the info of the supporting segment
    return info0;
  }
};

// functor that defines how to merge color info when a site (either
// point or segment) corresponds to point(s) on plane belonging to
// more than one input site
struct Red_blue_merge_info
{
  typedef Red_blue   Info;
  typedef Info       result_type;

  inline
  Info operator()(const Info& info0, const Info& info1) const {
    // if the two sites defining the new site have the same info, keep
    // this common info
    if ( info0 == info1 ) { return info0; }
    // otherwise the new site should be purple
    return -1;
  }
};

// define the storage traits with info
typedef CGAL::Segment_Delaunay_graph_storage_traits_with_info_2<Gt, Red_blue, Red_blue_convert_info, Red_blue_merge_info> ST;

typedef CGAL::Segment_Delaunay_graph_2<Gt,ST>  SDG;

#include <CGAL/IO/WKT.h>
#include <CGAL/Multipolygon_with_holes_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Segment_2.h>
using Point_2 = Kernel::Point_2;
using Polygon_2 = CGAL::Polygon_2<Kernel>;
using Segment_2 = CGAL::Segment_2<Kernel>;
using Multipolygon_with_holes_2 = CGAL::Multipolygon_with_holes_2<Kernel>;
using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<Kernel>;

using namespace std;

std::string currentDateTime() {
    std::time_t t = std::time(nullptr);
    std::tm* now = std::localtime(&t);
 
    char buffer[128];
    strftime(buffer, sizeof(buffer), "%m-%d-%Y %X", now);
    return buffer;
}
///////////////////////////////// CODE ABOUT EXACT DUALS ///////////////////////////////////////////

template < typename ExactSite, typename ExactKernel,
           typename SDGSite,   typename InputKernel >
ExactSite convert_site_to_exact(const SDGSite &site,
                                const InputKernel & /*k*/,
                                const ExactKernel & /*ek*/)
{
  using To_exact = CGAL::Cartesian_converter<InputKernel, ExactKernel>;
  To_exact to_exact;

  // Note: in theory, a site can be constructed from more than just one or two points
  // (e.g. 4 points for the segment defined by the intersection of two segments). Thus, it
  // would be better to convert the input points at the very beginning and just maintain
  // a type of map between the base and exact sites.
  ExactSite es;
  if(site.is_point())
    es = ExactSite::construct_site_2(to_exact(site.point()));
  else
    es = ExactSite::construct_site_2(to_exact(site.segment().source()), to_exact(site.segment().target()));

  return es;
}

// Dual (Voronoi site) of an SDG face
template < typename FiniteFacesIterator, typename InputKernel, typename ExactKernel >
typename ExactKernel::Point_2 exact_primal(const FiniteFacesIterator sdg_f,
                                           const InputKernel& k,
                                           const ExactKernel& ek)
{
  using Exact_SDG_traits = CGAL::Segment_Delaunay_graph_traits_2<ExactKernel>;
  using Exact_site_2 = typename Exact_SDG_traits::Site_2;

  static Exact_SDG_traits e_sdg_gt;
  const Exact_site_2 es0 = convert_site_to_exact<Exact_site_2>(sdg_f->vertex(0)->site(), k, ek);
  const Exact_site_2 es1 = convert_site_to_exact<Exact_site_2>(sdg_f->vertex(1)->site(), k, ek);
  const Exact_site_2 es2 = convert_site_to_exact<Exact_site_2>(sdg_f->vertex(2)->site(), k, ek);

  return e_sdg_gt.construct_svd_vertex_2_object()(es0, es1, es2);
}

// Dual (Voronoi edge) of an SDG edge
// this function is identical 'SDG::primal()', but with a conversion to exact sites
template < typename Edge, typename InputKernel, typename ExactKernel >
CGAL::Object exact_primal(const Edge& e,
                          const SDG& sdg,
                          const InputKernel& k,
                          const ExactKernel& ek)
{
  using Exact_SDG_traits = CGAL::Segment_Delaunay_graph_traits_2<ExactKernel>;
  using Exact_site_2 = typename Exact_SDG_traits::Site_2;

  using DT = CGAL::Field_with_sqrt_tag;
  using Construct_sdg_bisector_2 = CGAL::SegmentDelaunayGraph_2::Construct_sdg_bisector_2<Exact_SDG_traits, DT>;
  using Construct_sdg_bisector_ray_2 = CGAL::SegmentDelaunayGraph_2::Construct_sdg_bisector_ray_2<Exact_SDG_traits, DT>;
  using Construct_sdg_bisector_segment_2 = CGAL::SegmentDelaunayGraph_2::Construct_sdg_bisector_segment_2<Exact_SDG_traits, DT>;

  CGAL_precondition(!sdg.is_infinite(e));

  if(sdg.dimension() == 1)
  {
    Exact_site_2 p = convert_site_to_exact<Exact_site_2>((e.first)->vertex(sdg.cw(e.second))->site(), k, ek);
    Exact_site_2 q = convert_site_to_exact<Exact_site_2>((e.first)->vertex(sdg.ccw(e.second))->site(), k, ek);

    return make_object(Construct_sdg_bisector_2()(p, q));
  }

  // dimension == 2
  // neither of the two adjacent faces is infinite
  if((!sdg.is_infinite(e.first)) && (!sdg.is_infinite(e.first->neighbor(e.second))))
  {
    Exact_site_2 p = convert_site_to_exact<Exact_site_2>((e.first)->vertex(sdg.ccw(e.second))->site(), k, ek);
    Exact_site_2 q = convert_site_to_exact<Exact_site_2>((e.first)->vertex(sdg.cw(e.second))->site(), k, ek);
    Exact_site_2 r = convert_site_to_exact<Exact_site_2>((e.first)->vertex(e.second)->site(), k, ek);
    Exact_site_2 s = convert_site_to_exact<Exact_site_2>(sdg.tds().mirror_vertex(e.first, e.second)->site(), k, ek);

    return Construct_sdg_bisector_segment_2()(p, q, r, s);
  }

  // both of the adjacent faces are infinite
  if(sdg.is_infinite(e.first) && sdg.is_infinite(e.first->neighbor(e.second)))
  {
    Exact_site_2 p = convert_site_to_exact<Exact_site_2>((e.first)->vertex(sdg.cw(e.second))->site(), k, ek);
    Exact_site_2 q = convert_site_to_exact<Exact_site_2>((e.first)->vertex(sdg.ccw(e.second))->site(), k, ek);

    return make_object(Construct_sdg_bisector_2()(p, q));
  }

  // only one of the adjacent faces is infinite
  CGAL_assertion(sdg.is_infinite(e.first) || sdg.is_infinite(e.first->neighbor(e.second)));
  CGAL_assertion(!(sdg.is_infinite(e.first) && sdg.is_infinite(e.first->neighbor(e.second))));
  CGAL_assertion(sdg.is_infinite(e.first->vertex(e.second)) || sdg.is_infinite(sdg.tds().mirror_vertex(e.first, e.second)));

  Edge ee = e;
  if(sdg.is_infinite(e.first->vertex(e.second)))
  {
    ee = Edge(e.first->neighbor(e.second),
              e.first->neighbor(e.second)->index(sdg.tds().mirror_vertex(e.first, e.second)));
  }

  Exact_site_2 p = convert_site_to_exact<Exact_site_2>(ee.first->vertex(sdg.ccw(ee.second))->site(), k, ek);
  Exact_site_2 q = convert_site_to_exact<Exact_site_2>(ee.first->vertex(sdg.cw(ee.second))->site(), k, ek);
  Exact_site_2 r = convert_site_to_exact<Exact_site_2>(ee.first->vertex(ee.second)->site(), k, ek);

  return make_object(Construct_sdg_bisector_ray_2()(p, q, r));
}

///////////////////////////////// CODE TO DRAW A CROPPED DIAGRAM ///////////////////////////////////

template< typename OutputKernel>
void export_segment(const CGAL::Segment_2<OutputKernel>& segment, Red_blue info, const double xmin, const double ymin, std::ofstream& ofile)
{
  ofile << "\"LINESTRING(";
  ofile << CGAL::to_double(segment.source().x()+xmin) << " ";
  ofile << CGAL::to_double(segment.source().y()+ymin) << ", ";
  ofile << CGAL::to_double(segment.target().x()+xmin) << " ";
  ofile << CGAL::to_double(segment.target().y()+ymin) << ")\",";
  ofile << info << "\n" << std::flush;
}

// Split a Voronoi edge that is a parabola (one site is a point, one site is a segment) into small segments
template <typename OutputKernel>
void segment_parabola(const CGAL::Parabola_segment_2<OutputKernel>& p,
                      const CGAL::Bbox_2& scaled_bbox, Red_blue info, const double xmin, const double ymin, std::ofstream& ofile)
{
  using FT = typename OutputKernel::FT;
  using Point_2 = typename OutputKernel::Point_2;
  using Segment_2 = typename OutputKernel::Segment_2;

  const Point_2& o = p.origin();
  const Point_2& c = p.center();

  // @todo could be cached
  const Point_2 mm(scaled_bbox.xmin()+xmin, scaled_bbox.ymin()+ymin);
  const Point_2 mM(scaled_bbox.xmin()+xmin, scaled_bbox.ymax()+ymin);
  const Point_2 Mm(scaled_bbox.xmax()+xmin, scaled_bbox.ymin()+ymin);
  const Point_2 MM(scaled_bbox.xmax()+xmin, scaled_bbox.ymax()+ymin);

  FT s = CGAL::squared_distance(mm, c);
  s = (std::max)(s, CGAL::squared_distance(mM, c));
  s = (std::max)(s, CGAL::squared_distance(Mm, c));
  s = (std::max)(s, CGAL::squared_distance(MM, c));

  s = CGAL::sqrt(s) - CGAL::sqrt(CGAL::squared_distance(c, o));

  const double lx = scaled_bbox.xmax() - scaled_bbox.xmin();
  const double ly = scaled_bbox.ymax() - scaled_bbox.ymin();

  const double max_length = 10.0;//0.001 * (std::min)(lx, ly); // max length in the discretization of a parabola

  std::vector<Point_2> points;
  p.generate_points(points, max_length, -s, s);

  if(points.size() < 2) return;

  for(std::size_t i=0, ps=points.size()-1; i<ps; ++i)
    export_segment(Segment_2(points[i], points[i+1]), info, xmin, ymin, ofile);
    //segment_list.emplace_back(points[i], points[i+1]);
}

template< typename OutputKernel>
void fill_Voronoi_structure(const SDG& sdg,
                            const CGAL::Bbox_2& scaled_bbox,
                            const double xmin, const double ymin, 
                            std::ofstream& ofile)
{
  using Line_2 = typename OutputKernel::Line_2;
  using Ray_2 = typename OutputKernel::Ray_2;
  using Segment_2 = typename OutputKernel::Segment_2;
  using SDG_traits = CGAL::Segment_Delaunay_graph_traits_2<OutputKernel>;

  CK k;
  OutputKernel ek;

  Line_2 l;
  Ray_2 r;
  Segment_2 s;
  CGAL::Parabola_segment_2<SDG_traits> p;

  int nl = 0, ns = 0, nr = 0, np = 0;

  typename SDG::Finite_edges_iterator eit = sdg.finite_edges_begin(),
                                      eend = sdg.finite_edges_end();
  for (; eit != eend; ++eit)
  {
    CGAL::Object o = exact_primal(*eit, sdg, k, ek);

    Red_blue r_info = eit->first->vertex( sdg.ccw(eit->second) )->storage_site().info();
    Red_blue s_info = eit->first->vertex( sdg.cw(eit->second) )->storage_site().info();
    Red_blue info = (r_info == s_info) ? r_info : -1;

    if(CGAL::assign(l, o)) { /*line_list.push_back(l);*/ ++nl; }
    if(CGAL::assign(s, o)) { /*segment_list.push_back(s);*/ export_segment(s,info,xmin,ymin,ofile); ++ns; }
    if(CGAL::assign(r, o)) { /*ray_list.push_back(r);*/ ++nr; }
    if(CGAL::assign(p, o)) { segment_parabola(p, scaled_bbox, info, xmin,ymin, ofile); ++np; }
  }
}

Point_2 translate(const Point_2& p, double x, double y)
{
  return Point_2(p.x()+x, p.y()+y);
}

int main()
{
  cout << currentDateTime() << ": start" << endl;
  ifstream in("face.wkt");
  Multipolygon_with_holes_2 mp;
  CGAL::IO::read_multi_polygon_WKT(in, mp);
  cout << currentDateTime() << ": computing bbox" << endl;
  std::set<Point_2> all_points;
  for(Polygon_with_holes_2 p : mp) { 
    for(const Segment_2& e  : p.outer_boundary().edges()){
      all_points.insert(e.source());
      all_points.insert(e.target());
    }
  }

  // Get the bbox of the input points, and grow it a bit
  const CGAL::Bbox_2 bbox = bbox_2(all_points.begin(), all_points.end());
  const double xmin = bbox.xmin(), xmax = bbox.xmax();
  const double ymin = bbox.ymin(), ymax = bbox.ymax();
  const double xmid = 0.5 * (xmin + xmax), ymid = 0.5 * (ymin + ymax);
  const double scaling_factor = 3.; // '0.5' gives the identity
  const double lx = scaling_factor * (xmax - xmin),
               ly = scaling_factor * (ymax - ymin);
  const CGAL::Bbox_2 scaled_bbox(xmid - lx, ymid - ly, xmid + lx, ymid + ly);

  std::cout << "bbox: " << bbox.xmin() << " " << bbox.ymin() << std::endl;
  std::cout << "bbox: " << bbox.xmax() << " " << bbox.ymax() << std::endl;
  std::cout << "lx/y: " << lx << " " << ly << std::endl;
  std::cout << "Scaled bbox: " << scaled_bbox.xmin() << " " << scaled_bbox.ymin() << std::endl;
  std::cout << "Scaled bbox: " << scaled_bbox.xmax() << " " << scaled_bbox.ymax() << std::endl;

  SDG          sdg;
  SDG::Site_2  site;
  vector<SDG::Site_2> sites;

  int ring_index = 0;
  for(Polygon_with_holes_2 p : mp) { 
    cout << "polygon " <<ring_index<< endl;
    for(const Segment_2& e  : p.outer_boundary().edges()){
      //cout << e << endl;
      site = SDG::Site_2::construct_site_2(translate(e.source(),-xmin,-ymin), translate(e.target(),-xmin,-ymin));
      //site.set_info(ring_index);
      //sites.push_back(site);
      sdg.insert(site, ring_index);
    }
    ring_index++;
    typename Polygon_with_holes_2::Hole_const_iterator  hit;
    for (hit = p.holes_begin(); hit != p.holes_end(); ++hit) {
      cout << "hole " <<ring_index<< endl;
      for(const Segment_2& e  : hit->edges()){
        //cout << e << endl;
        site = SDG::Site_2::construct_site_2(translate(e.source(),-xmin,-ymin), translate(e.target(),-xmin,-ymin));
        //site.set_info(ring_index);
        //sites.push_back(site);
        sdg.insert(site, ring_index);
      }
      ring_index++;
    }
  }
  cout << currentDateTime() << ": inserting " << ring_index << " rings" <<endl;
  //insert the sites all at once using spatial sorting to speed the insertion
  //sdg.insert( sites.begin(), sites.end(), CGAL::Tag_true() );
 
  // validate the diagram
  //assert( sdg.is_valid(true, 1) );
  cout << sdg.is_valid(true, 1) << endl;
  // Output to WKT file
  std::ofstream contour_ofile ("sdg_info.wkt");
  contour_ofile.precision(18);
  contour_ofile << "wkt,ring\n";
  cout << currentDateTime() << ": fill_Voronoi_structure" << endl;
  fill_Voronoi_structure<EK>(sdg, scaled_bbox, xmin, ymin, contour_ofile);
  contour_ofile.close();
  cout << currentDateTime() << ": done" << endl;
  return 0;
}
