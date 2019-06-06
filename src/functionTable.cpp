#include <astro/utilities/functionTable.h>
#include <astro/utilities/utilities.h>

astro::functionTable::functionTable ():
nx (0), ny (0), cut_min (false), cut_max (false)
{
  xt.clear ();
  yt.clear ();
  vt.clear ();
  at.clear ();
}

astro::functionTable::functionTable
  (const std::vector<double>& x, const std::vector<double>& v):
nx (x.size ()), ny (0), cut_min (false), cut_max (false), xt (x), vt (v)
{
  if (nx != int (vt.size ()))
    throw std::out_of_range
      ("in functionTable::fill: size mismatch in input arrays");
  yt.clear ();
  at.clear ();
}

astro::functionTable::functionTable
  (const std::vector<double>& x, const std::vector<double>& y,
   const std::vector< std::vector<double> >& a):
nx (x.size ()), ny (y.size ()), cut_min (false), cut_max (false),
xt (x), yt (y), at (a)
{
  if (nx != int (at.size ()) || ny != int (at[0].size ()))
    throw std::out_of_range
      ("in functionTable::fill: size mismatch in input arrays");
  vt.clear ();
}

astro::functionTable::functionTable
  (std::function<double (double)> fct, const astro::spacing_type type,
   const unsigned int n, const double min, const double max):
nx (n), ny (0), cut_min (false), cut_max (false)
{
  xt = astro::fill (type, nx, min, max);
  for (int i = 0; i < nx; i++)
    vt.push_back (fct (xt[i]));
  at.clear ();
}

astro::functionTable::functionTable
  (std::function<double (double, double)> fct,
   const spacing_type x_type, const spacing_type y_type,
   const unsigned int nx_in, const unsigned int ny_in,
   const double x_min, const double x_max,
   const double y_min, const double y_max):
nx (nx_in), ny (ny_in), cut_min (false), cut_max (false)
{
  xt = astro::fill (x_type, nx, x_min, x_max);
  yt = astro::fill (y_type, ny, y_min, y_max);
  vt.clear ();
  at.resize (nx);
  for (int i = 0; i < nx; i++)
  {
    at[i].resize (ny);
    for (int j = 0; j < ny; j++)
      at[i][j] = fct (xt[i], yt[j]);
  }
}

void astro::functionTable::set_cut_min (bool cut_min_in, double null_min_in)
{
  cut_min = cut_min_in;
  null_min = null_min_in;
}

void astro::functionTable::set_cut_max (bool cut_max_in, double null_max_in)
{
  cut_max = cut_max_in;
  null_max = null_max_in;
}

void astro::functionTable::set_cut
  (bool cut_in, double null_min_in, double null_max_in)
{
  cut_min = cut_in;
  cut_max = cut_in;
  null_min = null_min_in;
  null_max = null_max_in;
}

double astro::functionTable::operator() (const double x) const
{
  if (ny > 0)
    throw std::invalid_argument
      ("in functionTable::operator (): function needs two arguments");
  if (cut_min && x < xt[0])
    return null_min;
  if (cut_max && x > xt[nx-1])
    return null_max;
  int ix = astro::locate (xt, x);
  ix = std::min (std::max (ix,0), nx-2);
  double fx = (x-xt[ix])/(xt[ix+1]-xt[ix]);
  return
     fx     *vt[ix+1]+
    (1.0-fx)*vt[ix];
}

double astro::functionTable::operator () (const double x, const double y) const
{
  if (ny == 0)
    throw std::invalid_argument
      ("in functionTable::operator (): function needs one argument");
  if (cut_min && (x < xt[0] || y < yt[0]))
    return null_min;
  if (cut_max && (x > xt[nx-1] || y > yt[ny-1]))
    return null_max;
  int ix = astro::locate (xt, x);
  ix = std::min (std::max (ix,0), nx-2);
  int iy = astro::locate (yt, y);
  iy = std::min (std::max (iy,0), ny-2);
  double fx = (x-xt[ix])/(xt[ix+1]-xt[ix]);
  double fy = (y-yt[iy])/(yt[iy+1]-yt[iy]);
  return
     fx     *     fy *at[ix+1][iy+1]+
    (1.0-fx)*     fy *at[ix  ][iy+1]+
     fx     *(1.0-fy)*at[ix+1][iy  ]+
    (1.0-fx)*(1.0-fy)*at[ix  ][iy  ];
}

void astro::functionTable::fill
  (const std::vector<double>& x, const std::vector<double>& v)
{
  if (x.size () != v.size ())
    throw std::out_of_range
      ("in functionTable::fill: size mismatch in input arrays");
  nx = x.size ();
  ny = 0;
  xt = x;
  vt = v;
  at.clear ();
}

void astro::functionTable::fill
  (const std::vector<double>& x, const std::vector<double>& y,
   const std::vector< std::vector<double> >& a)
{
  if (x.size () != a.size () || y.size () != a[0].size ())
    throw std::out_of_range
      ("in functionTable::fill: size mismatch in input arrays");
  nx = x.size ();
  ny = y.size ();
  xt = x;
  yt = y;
  vt.clear ();
  at = a;
}
