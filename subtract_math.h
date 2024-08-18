#if !defined(SUBTRACT_MATH_H)
/* ========================================================================
   $File: $
   $Date: 2024 $
   $Revision: $
   $Creator: BabyKaban $
   $Notice: $
   ======================================================================== */
#include "math.h"

inline v2
V2(f32 X, f32 Y)
{
    v2 Result;

    Result.x = X;
    Result.y = Y;

    return(Result);
}

inline f32
SquareRoot(f32 F32)
{
    f32 Result = sqrtf(F32);
    return(Result);
}

inline f32
AbsoluteValue(f32 F32)
{
    f32 Result = (f32)fabs(F32);
    return(Result);
}

inline f32
Lerp(f32 A, f32 t, f32 B)
{
    f32 Result = (1.0f - t)*A + t*B;

    return(Result);
}

inline f32
Clamp(f32 Min, f32 Value, f32 Max)
{
    f32 Result = Value;

    if(Result < Min)
    {
        Result = Min;
    }
    else if(Result > Max)
    {
        Result = Max;
    }

    return(Result);
}

inline f32
Clamp01(f32 Value)
{
    f32 Result = Clamp(0.0f, Value, 1.0f);

    return(Result);
}

inline v2
Perp(v2 A)
{
    v2 Result = {-A.y, A.x};
    return(Result);
}

inline v2
operator*(f32 A, v2 B)
{
    v2 Result;

    Result.x = A*B.x;
    Result.y = A*B.y;
    
    return(Result);
}

inline v2
operator*(v2 B, f32 A)
{
    v2 Result = A*B;

    return(Result);
}

inline v2 &
operator*=(v2 &B, f32 A)
{
    B = A * B;

    return(B);
}

inline v2
operator-(v2 A)
{
    v2 Result;

    Result.x = -A.x;
    Result.y = -A.y;

    return(Result);
}

inline v2
operator+(v2 A, v2 B)
{
    v2 Result;

    Result.x = A.x + B.x;
    Result.y = A.y + B.y;

    return(Result);
}

inline v2 &
operator+=(v2 &A, v2 B)
{
    A = A + B;

    return(A);
}

inline v2
operator-(v2 A, v2 B)
{
    v2 Result;

    Result.x = A.x - B.x;
    Result.y = A.y - B.y;

    return(Result);
}

inline v2 &
operator-=(v2 &A, v2 B)
{
    A = A - B;

    return(A);
}

inline v2
Lerp(v2 A, f32 t, v2 B)
{
    v2 Result = (1.0f - t)*A + t*B;

    return(Result);
}

inline v2
Hadamard(v2 A, v2 B)
{
    v2 Result = {A.x*B.x, A.y*B.y};

    return(Result);
}

inline f32
Inner(v2 A, v2 B)
{
    f32 Result = A.x*B.x + A.y*B.y;

    return(Result);
}

inline f32
Cross(v2 A, v2 B)
{
    f32 Result = A.x*B.y - A.y*B.x;

    return(Result);
}

inline f32
LengthSq(v2 A)
{
    f32 Result = Inner(A, A);

    return(Result);
}

inline f32
Length(v2 A)
{
    f32 Result = SquareRoot(LengthSq(A));
    return(Result);
}

inline v2
Clamp01(v2 Value)
{
    v2 Result;

    Result.x = Clamp01(Value.x);
    Result.y = Clamp01(Value.y);

    return(Result);
}

inline v2
Normalize(v2 A)
{
    v2 Result = {};

    f32 L = Length(A);
    if(L)
    {
        Result = A * (1.0f / L);
    }

    return(Result);
}

struct polygon2
{
    s32 VertexCount;
    v2 *Vertices;

    b32 HasHole;
    s32 HoleVertexCount;
    v2 *HoleVertices;
};

#define MAX_POLYGON_COUNT 32
struct polygon2_set
{
    s32 PolygonCount;
    polygon2 *Polygons;
};

struct triangle
{
    v2 Vertices[3];
};

inline f32
TriangleSignedArea(v2 a, v2 b, v2 c)
{
    f32 Result = 0.5f*(a.x*(b.y - c.y) + b.x*(c.y - a.y) + c.x*(a.y - b.y));
    return(Result);
}

inline f32
TriangleSignedArea(triangle *T)
{
    f32 Result = TriangleSignedArea(T->Vertices[0], T->Vertices[1], T->Vertices[2]);
    return(Result);
}

inline f32
PolygonSignedArea(polygon2 *Polygon)
{
    f32 Result = 0.0f;
    for(s32 Index = 0;
        Index < Polygon->VertexCount;
        ++Index)
    {
        s32 Next = (Index + 1) % Polygon->VertexCount;
        v2 FirstVertex = Polygon->Vertices[Index];
        v2 SecondVertex = Polygon->Vertices[Next];

        Result += Cross(FirstVertex, SecondVertex);
    }

    Result = 0.5f*Result;
    
    return(Result);
}

#define SUBTRACT_MATH_H
#endif
