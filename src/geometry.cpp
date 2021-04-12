#include "geometry.h"

template <> template <> vec<3,int>  ::vec(const vec<3,double> &v) : x(int(v.x+.5)),y(int(v.y+.5)),z(int(v.z+.5)) {}
template <> template <> vec<3,double>::vec(const vec<3,int> &v)   : x(v.x),y(v.y),z(v.z) {}
template <> template <> vec<2,int>  ::vec(const vec<2,double> &v) : x(int(v.x+.5)),y(int(v.y+.5)) {}
template <> template <> vec<2,double>::vec(const vec<2,int> &v)   : x(v.x),y(v.y) {}