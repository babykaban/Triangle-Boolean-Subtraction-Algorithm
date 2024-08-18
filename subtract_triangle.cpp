/* ========================================================================
   $File: $
   $Date: 2024 $
   $Revision: $
   $Creator: BabyKaban $
   $Notice:  $
   ======================================================================== */
#include <stdio.h>
#include <stdlib.h>

#include "subtract_triangle.h"
#include "subtract_math.h"

#define SUBTRACTION_MAX_POINTS_PER_POLYGON 16 // NOTE(babykaban): Actually 10 but to make sure that there is enough spece 16 should be good
#define SUBTRACTION_MAX_POLYGON_COUNT 4 // NOTE(babykaban): Actually 3 
#define SUBTRACT_EPSILON_F32 0.000001f
#define SUBTRACT_MINIMAL_POINT_DISTANCE 0.00008f

struct subtract_result
{
    polygon2_set Set;
    b32 Success;
};

inline void
ProjectVerticesOntoNormal(v2 *Vertices, v2 Normal, f32 *Min, f32 *Max)
{
    *Min = *Max = Inner(Vertices[0], Normal);
    for(s32 I = 1;
        I < 3;
        ++I)
    {
        f32 Projection = Inner(Vertices[I], Normal);
        if(Projection < *Min) *Min = Projection;
        if(Projection > *Max) *Max = Projection;
    }
}

inline b32
OverlapExcludingVerticesAndEdges(f32 Min1, f32 Max1, f32 Min2, f32 Max2)
{
    b32 Result = !((Max1 < Min2) || (Max2 < Min1) ||
                   (Max1 == Min2) || (Max2 == Min1));

    return(Result);
}

internal b32
TrianglesOverlapExcludingVerticesAndEdges(triangle *A, triangle *B)
{
    b32 Result = true;

    v2 Edges1[3] = 
        {
            A->Vertices[1] - A->Vertices[0],
            A->Vertices[2] - A->Vertices[1],
            A->Vertices[0] - A->Vertices[2],
        };

    v2 Edges2[3] =
        {
            B->Vertices[1] - B->Vertices[0],
            B->Vertices[2] - B->Vertices[1],
            B->Vertices[0] - B->Vertices[2],
        };

    v2 Normals[6];
    for(s32 I = 0;
        I < 3;
       ++I)
    {
        Normals[I] = Perp(Edges1[I]);
        Normals[I + 3] = Perp(Edges2[I]);
    }

    for(s32 I = 0;
        I < 6;
        ++I)
    {
        f32 Min1, Max1, Min2, Max2;

        ProjectVerticesOntoNormal(A->Vertices, Normals[I], &Min1, &Max1);
        ProjectVerticesOntoNormal(B->Vertices, Normals[I], &Min2, &Max2);

        if(!OverlapExcludingVerticesAndEdges(Min1, Max1, Min2, Max2))
        {
            Result = false;
            break;
        }
    }

    return(Result);
}

inline b32
PointsAreEqualF32(v2 a, v2 b)
{
    b32 Result = ((AbsoluteValue(a.x - b.x) < SUBTRACT_EPSILON_F32) && (AbsoluteValue(a.y - b.y) < SUBTRACT_EPSILON_F32));
    return(Result);
}

inline b32
IsNewPoint(polygon2 *A, v2 p)
{
    b32 Result = true;
    for(s32 I = 0;
        I < A->VertexCount;
        ++I)
    {
        if(PointsAreEqualF32(A->Vertices[I], p))
        {
            Result = false;
            break;
        }
    }

    return(Result);
}

inline b32
IsNewPoint(v2 *A, s32 Count, v2 p)
{
    b32 Result = true;
    for(s32 I = 0;
        I < Count;
        ++I)
    {
        if(PointsAreEqualF32(A[I], p))
        {
            Result = false;
            break;
        }
    }

    return(Result);
}

inline b32
LineSectF32(v2 x0, v2 x1, v2 y0, v2 y1, v2 *res)
{
    b32 Result = true;

    v2 dx = x1 - x0;
    v2 dy = y1 - y0;
    v2 d = x0 - y0;

    f32 dyx = Cross(dy, dx);
    if(AbsoluteValue(dyx) < SUBTRACT_EPSILON_F32)
    {
        Result = false;
    }
    else
    {
        dyx = Cross(d, dx) / dyx;
        if(dyx <= 0 || dyx >= 1)
        {
            Result = false;
        }
        else
        {
            *res = y0 + dyx * dy;
        }
    }

    return(Result);
}

inline f32
DistanceToSegment(v2 p, v2 a, v2 b)
{
    f32 Result = 0.0f;
    
    f32 l2 = LengthSq(a - b);
    if(l2 == 0.0f)
    {
        Result = Length(p - a);
    }
    else
    {
        v2 pa = p - a;
        v2 ba = b - a;
        f32 t = Inner(pa, ba) / l2;
        t = Clamp(0, t, 1);

        v2 Closest = Lerp(a, t, b);

        Result = Length(p - Closest);
    }

    return(Result);
}

inline b32
IsTriangleCollinear(triangle *A, f32 Epsilon)
{
    f32 d0 = DistanceToSegment(A->Vertices[0], A->Vertices[1], A->Vertices[2]);
    f32 d1 = DistanceToSegment(A->Vertices[1], A->Vertices[0], A->Vertices[2]);
    f32 d2 = DistanceToSegment(A->Vertices[2], A->Vertices[0], A->Vertices[1]);

    b32 Result = ((d0 < Epsilon) || (d1 < Epsilon) || (d2 < Epsilon));
    return(Result);
}

inline b32
InsideF32(v2 p, v2 Clipper1, v2 Clipper2)
{
    b32 Result = false;
    f32 CrossProduct = (Clipper2.x - Clipper1.x) * (p.y - Clipper1.y) - (Clipper2.y - Clipper1.y) * (p.x - Clipper1.x);
    Result = CrossProduct < SUBTRACT_EPSILON_F32;

    return(Result);
}

inline s32
LeftOfF32(v2 a, v2 b, v2 c)
{
    s32 Result = 0;

    v2 ba = b - a;
    v2 cb = c - b;

    f32 x = Cross(ba, cb);

    Result = (x < SUBTRACT_EPSILON_F32) ? -1 : (x > 0);
    return(Result);
}

inline void
PolyEdgeClipF32(polygon2 *Sub, v2 x0, v2 x1, s32 Left, polygon2 *Res)
{
    v2 v0 = Sub->Vertices[Sub->VertexCount- 1];
    v2 v1 = {};

    Res->VertexCount = 0;

    s32 Side0 = LeftOfF32(x0, x1, v0);
    if(Side0 != -Left)
    {
        Res->Vertices[Res->VertexCount++] = v0;
    }

    for(s32 I = 0;
        I < Sub->VertexCount;
        ++I)
    {
        v1 = Sub->Vertices[I];
        s32 Side1 = LeftOfF32(x0, x1, v1);
        if(((Side0 + Side1) == 0) && Side0)
        {
            v2 Intersect;
            if(LineSectF32(x0, x1, v0, v1, &Intersect))
            {
                Res->Vertices[Res->VertexCount++] = Intersect;
            }
        }

        if(I == (Sub->VertexCount - 1))
        {
            break;
        }

        if(Side1 != -Left)
        {
            Res->Vertices[Res->VertexCount++] = v1;
        }

        v0 = v1;
        Side0 = Side1;
    }
}

internal void
SutherlandHodgman(polygon2 *Subject, triangle *Clip, polygon2 *Result)
{
    polygon2 Temp;
    Temp.VertexCount = 0;
    Temp.Vertices = (v2 *)malloc(SUBTRACTION_MAX_POINTS_PER_POLYGON*sizeof(v2));

    polygon2 New0;
    New0.VertexCount = 0;
    New0.Vertices = (v2 *)malloc(SUBTRACTION_MAX_POINTS_PER_POLYGON*sizeof(v2));

    s32 dir = LeftOfF32(Clip->Vertices[0], Clip->Vertices[1], Clip->Vertices[2]);
    PolyEdgeClipF32(Subject, Clip->Vertices[2], Clip->Vertices[0], dir, Result);
    for(s32 I = 0;
        I < 2;
        ++I)
    {
        Temp.VertexCount = Result->VertexCount;
        Copy(sizeof(v2)*Result->VertexCount, Result->Vertices, Temp.Vertices);

        Result->VertexCount = New0.VertexCount;
        Copy(sizeof(v2)*New0.VertexCount, New0.Vertices, Result->Vertices);

        New0.VertexCount = Temp.VertexCount;
        Copy(sizeof(v2)*Temp.VertexCount, Temp.Vertices, New0.Vertices);

        if(Temp.VertexCount == 0)
        {
            Result->VertexCount = 0;
            break;
        }

        PolyEdgeClipF32(&New0, Clip->Vertices[I], Clip->Vertices[I + 1], dir, Result);
    }

    free(Temp.Vertices);
    free(New0.Vertices);
}

inline void
CleanUpOverlapPolygon(polygon2 *Poly)
{
    for(s32 I = 0;
        I < Poly->VertexCount;
        ++I)
    {
        v2 Cur = Poly->Vertices[I];
        v2 Next = Poly->Vertices[(I + 1) % Poly->VertexCount];
        f32 l = Length(Cur - Next);
        if(l < SUBTRACT_MINIMAL_POINT_DISTANCE)
        {
            Poly->Vertices[I + 1] = {};
            for(s32 J = I + 1;
                J < (Poly->VertexCount - 1);
                ++J)
            {
                Poly->Vertices[J] = Poly->Vertices[J + 1];                
            }
            --Poly->VertexCount;
            --I;
        }
    }

    Poly->Vertices[Poly->VertexCount] = {};
}

internal void
HandleVertexCase(polygon2 *A, polygon2 *B)
{
    b32 RotateVerticies = false;
    for(s32 VertexIndex = 0;
        VertexIndex < B->VertexCount;
        ++VertexIndex)
    {
        if(PointsAreEqualF32(A->Vertices[0], B->Vertices[VertexIndex]))
        {
            RotateVerticies = true;
            break;
        }
    }
    
    if(RotateVerticies)
    {
        v2 Temp = A->Vertices[2];
        A->Vertices[2] = A->Vertices[1];
        A->Vertices[1] = A->Vertices[0];
        A->Vertices[0] = Temp;
    }
}

inline void
FindOutsidePointsForF32(polygon2 *A, polygon2 *B, v2 *OutsidePoints, s32 *OutsideCount)
{
    for(s32 VertextIndex = 0;
        VertextIndex < A->VertexCount;
        ++VertextIndex)
    {
        b32 Inside = false;
        v2 p = A->Vertices[VertextIndex];
        for(s32 VertextIndexB = 0;
            VertextIndexB < B->VertexCount;
            ++VertextIndexB)
        {
            v2 p0 = B->Vertices[VertextIndexB];
            if(PointsAreEqualF32(p, p0))
            {
                Inside = true;
                break;
            }
        }

        if(!Inside)
        {
            OutsidePoints[*OutsideCount] = p;
            (*OutsideCount)++;
        }
    }
}

int
Compare(const void *a, const void *b)
{
    sort_entry *A = (sort_entry *)a; 
    sort_entry *B = (sort_entry *)b; 
    if(A->SortKey < B->SortKey)
    {
        return(-1);
    }
    else if(A->SortKey > B->SortKey)
    {
        return(1);
    }
    else
    {
        return(0);
    }
}

internal b32
IntersectLineSegment(polygon2 *A, v2 p1, v2 p2, v2 *Points, s32 *PointCount)
{
    b32 Intersect = false;
    s32 TempPointCount = 0;
    v2 *TempPoints = (v2 *)malloc(SUBTRACTION_MAX_POINTS_PER_POLYGON*sizeof(v2));

    for(s32 I = 0;
        I < (A->VertexCount - 1);
        ++I)
    {
        s32 J = (I + 1) % (A->VertexCount - 1);
        v2 p3 = A->Vertices[I];
        v2 p4 = A->Vertices[J];

        b32 Ends = false;
        if(DistanceToSegment(p3, p1, p2) < SUBTRACT_EPSILON_F32)
        {
            TempPoints[TempPointCount++] = p3;
            Intersect = true;
            Ends = true;
        }

        if(DistanceToSegment(p4, p1, p2) < SUBTRACT_EPSILON_F32)
        {
            TempPoints[TempPointCount++] = p4;
            Intersect = true;
            Ends = true;
        }
    }

    if(TempPointCount > 0)
    {
        ZeroArray((*PointCount)*sizeof(v2), Points);
        if(TempPointCount == 1)
        {
            Points[(*PointCount)++] = TempPoints[0];
        }
        else
        {
            s32 SortCount = TempPointCount;
            sort_entry *SortArray = (sort_entry *)malloc(TempPointCount*sizeof(sort_entry));
            for(s32 I = 0;
                I < TempPointCount;
                ++I)
            {
                sort_entry *Entry = SortArray + I;
                Entry->Index = I;
                Entry->SortKey = Length(p1 - TempPoints[I]);
            }

            qsort(SortArray, SortCount, sizeof(sort_entry), Compare);

            v2 Last = TempPoints[SortArray[0].Index];
            Points[(*PointCount)++] = Last;
            for(s32 I = 1;
                I < TempPointCount;
                ++I)
            {
                v2 This = TempPoints[SortArray[I].Index];
                if(!PointsAreEqualF32(This, Last))
                {
                    Points[(*PointCount)++] = This;
                    Last = This;
                }
            }
        }
    }

    free(TempPoints);
    
    return(Intersect);
}

struct points_between
{
    s32 Count;
    s32 I0;
    s32 I1;
};

internal points_between *
GetIntersectionPointsF32(polygon2 *A, polygon2 *B, v2 *OutputPoints, s32 *OutputCount)
{
    points_between *Result = (points_between *)malloc(3*sizeof(points_between));

    s32 TempCount = 0;
    v2 *Temp = (v2 *)malloc(SUBTRACTION_MAX_POINTS_PER_POLYGON*sizeof(v2));

    for(s32 I = 0;
        I < A->VertexCount;
        ++I)
    {
        v2 Prv = A->Vertices[I];
        v2 Cur = A->Vertices[(I + 1) % A->VertexCount];
        if(!PointsAreEqualF32(Prv, Cur))
        {
            IntersectLineSegment(B, Prv, Cur, Temp, &TempCount);
            s32 AddedCount = 0;
            for(s32 J = 0;
                J < TempCount;
                ++J)
            {
                if(IsNewPoint(OutputPoints, *OutputCount, Temp[J]))
                {
                    OutputPoints[(*OutputCount)++] = Temp[J];
                    ++AddedCount;
                }
            }

            Result[I].Count = AddedCount;
            Result[I].I0 = I;
            Result[I].I1 = (I + 1) % A->VertexCount;
            
            TempCount = 0;
        }
    }

    free(Temp);
    
    return(Result);
}

struct vertex_info
{
    s32 VertexIndex;
    v2 Point;

    b32 Outside;
    s32 Cross;
    b32 Processed;
};

internal void
BuildVertexTableForF32(polygon2 *A, polygon2 *B, v2 *OutsidePoints, s32 OutsideCount, vertex_info *Table)
{
    for(s32 VertexIndex = 0;
        VertexIndex < A->VertexCount;
        ++VertexIndex)
    {
        v2 p = A->Vertices[VertexIndex];

        vertex_info *Info = Table + VertexIndex;
        *Info = {};
        Info->VertexIndex = VertexIndex;
        Info->Point = p;
        Info->Processed = false;

        for(s32 OutsideIndex = 0;
            OutsideIndex < OutsideCount;
            ++OutsideIndex)
        {
            v2 OutsidePoint = OutsidePoints[OutsideIndex];
            if(PointsAreEqualF32(p, OutsidePoint))
            {
                Info->Outside = true;
                Info->Cross = -1;
                break;
            }
        }

        if(!Info->Outside)
        {
            Info->Cross = -1;
        }
    }
}

internal void
RecordCrossingPointsF32(vertex_info *A, s32 CountA, vertex_info *B, s32 CountB)
{
    for(s32 I = 0;
        I < CountA;
        ++I)
    {
        vertex_info *VertexA = A + I;
        if(!VertexA->Outside)
        {
            for(s32 J = 0;
                J < CountB;
                ++J)
            {
                vertex_info *VertexB = B + J;
                if(PointsAreEqualF32(VertexA->Point, VertexB->Point))
                {
                    VertexA->Cross = J;
                    VertexB->Cross = I;
                }
            }
        }
    }
}

struct between_indecies
{
    s32 Index0;
    s32 Index1;
};

inline void
InsertPointBetween(polygon2 *A, v2 p, between_indecies Indecies)
{
    for(s32 I = A->VertexCount;
        I > Indecies.Index1;
        --I)
    {
        A->Vertices[I] = A->Vertices[I - 1];
    }

    ++A->VertexCount;
    A->Vertices[Indecies.Index1] = p;
}

internal b32
ConstructPolygonsFor(polygon2 *A, polygon2 *B, vertex_info *InfoA, vertex_info *InfoB, polygon2_set *Set)
{
    b32 Result = true;
    for(;;)
    {
        if(!Result)
        {
            break;
        }
        
        s32 IndexA = -1;
        for(s32 VertexInfoIndexA = 0;
            VertexInfoIndexA < A->VertexCount;
            ++VertexInfoIndexA)
        {
            vertex_info *PolygonAInfo = InfoA + VertexInfoIndexA;
            if(PolygonAInfo->Outside && !PolygonAInfo->Processed)
            {
                IndexA = PolygonAInfo->VertexIndex;
                polygon2 *Polygon = Set->Polygons + Set->PolygonCount;
                Polygon->VertexCount = 0;
                Polygon->Vertices = (v2 *)malloc(sizeof(v2)*SUBTRACTION_MAX_POINTS_PER_POLYGON);

                ++Set->PolygonCount;
                break;
            }
        }

        polygon2 *CurrentPolygon = Set->Polygons + (Set->PolygonCount - 1);
        if(IndexA != -1)
        {
            for(;;)
            {
                if(CurrentPolygon->VertexCount >= SUBTRACTION_MAX_POINTS_PER_POLYGON)
                {
                    Result = false;
                    break;
                }
                
                CurrentPolygon->Vertices[CurrentPolygon->VertexCount++] = InfoA[IndexA].Point; 
                InfoA[IndexA].Processed = true;

                // NOTE(babykaban): Check if polygon completed
                v2 FirstOutputPoint = v2(CurrentPolygon->Vertices[0]);
                v2 ThisPoint = InfoA[IndexA].Point;
                if(CurrentPolygon->VertexCount > 1)
                {
                    if(PointsAreEqualF32(FirstOutputPoint, ThisPoint))
                    {
                        break;
                    }
                }

                s32 Cross = InfoA[IndexA].Cross;
                if(Cross != -1)
                {
                    s32 Next = ((Cross - 1) + B->VertexCount) % B->VertexCount; 
                    if(InfoB[Next].Outside)
                    {
                        Cross = -1;
                    }
                    else
                    {
                        s32 CrossBack = InfoB[Next].Cross;
                        if(CrossBack == (IndexA + 1) % A->VertexCount)
                        {
                            Cross = -1;
                        }
                    }
                }

                if(Cross != -1)
                {
                    s32 IndexB = Cross;
                    for(;;)
                    {
                        if(CurrentPolygon->VertexCount >= SUBTRACTION_MAX_POINTS_PER_POLYGON)
                        {
                            Result = false;
                            break;
                        }

                        IndexB = ((IndexB - 1) + B->VertexCount) % B->VertexCount;
                        CurrentPolygon->Vertices[CurrentPolygon->VertexCount++] = InfoB[IndexB].Point; 

                        Cross = InfoB[IndexB].Cross;
                        if(Cross != -1)
                        {
                            s32 Next = (Cross + 1) % A->VertexCount;
                            if(!InfoA[Next].Outside)
                            {
                                Cross = -1;
                            }
                            else
                            {
                                s32 CrossBack = InfoA[Next].Cross;
                                if(CrossBack == ((IndexB - 1) + B->VertexCount) % B->VertexCount)
                                {
                                    Cross = -1;
                                }
                            }
                        }

                        if(Cross != -1)
                        {
                            IndexA = Cross + 1;
                            break;
                        }
                    }
                }
                else
                {
                    ++IndexA;
                }
            }
        }
        else
        {
            break;
        }
    }

    return(Result);
}
 
inline void
RemoveDublicatPoints(polygon2 *Poly)
{
    for(s32 I = 0;
        I < Poly->VertexCount;
        ++I)
    {
        v2 Cur = Poly->Vertices[I];
        v2 Next = Poly->Vertices[(I + 1) % Poly->VertexCount];
        if(PointsAreEqualF32(Cur, Next))
        {
            Poly->Vertices[I + 1] = {};
            for(s32 J = I + 1;
                J < (Poly->VertexCount - 1);
                ++J)
            {
                Poly->Vertices[J] = Poly->Vertices[J + 1];                
            }
            --Poly->VertexCount;
            --I;
        }
    }

    Poly->Vertices[Poly->VertexCount] = {};
}

inline b32
IsPolygonCollinearF32(polygon2 *A, f32 Epsilon)
{
    b32 Result = false;

    for(s32 I = 0;
        I < A->VertexCount;
        ++I)
    {
        s32 First = I;
        s32 Second = (I + 1) % A->VertexCount;
        s32 Third = (I + 2) % A->VertexCount;
        f32 d = DistanceToSegment(A->Vertices[First], A->Vertices[Second], A->Vertices[Third]);
        if(d < Epsilon)
        {
            ZeroArray(sizeof(v2)*A->VertexCount, A->Vertices);
            A->VertexCount = 0;
            Result = true;
            break;
        }
    }
    return(Result);
}

internal subtract_result
SubtractTriangels(triangle *Subject, triangle *Subtractor, f32 MinimalOverlapArea, f32 CollinearValue, f32 MinimalResultArea)
{
    subtract_result Result = {};
    Result.Set.PolygonCount = 0;
    Result.Set.Polygons = (polygon2 *)malloc(sizeof(polygon2)*SUBTRACTION_MAX_POLYGON_COUNT);

    b32 ExcludingOverlap = TrianglesOverlapExcludingVerticesAndEdges(Subject, Subtractor);
    f32 AreaA = TriangleSignedArea(Subject);
    f32 AreaB = TriangleSignedArea(Subtractor);

    // NOTE(babykaban): Changing orrientation to be clockwise
    if(AreaA > 0.0f)
    {
        v2 Temp = Subject->Vertices[1];
        Subject->Vertices[1] = Subject->Vertices[2];
        Subject->Vertices[2] = Temp;
    }

    if(AreaB > 0.0f)
    {
        v2 Temp = Subtractor->Vertices[1];
        Subtractor->Vertices[1] = Subtractor->Vertices[2];
        Subtractor->Vertices[2] = Temp;
    }

    // NOTE(babykaban): Checks to ensure that both triangles are valid
    if(ExcludingOverlap && !IsTriangleCollinear(Subject, CollinearValue) && !IsTriangleCollinear(Subtractor, CollinearValue))
    {
        polygon2 PolygonA = {};
        PolygonA.Vertices = (v2 *)malloc(SUBTRACTION_MAX_POINTS_PER_POLYGON*sizeof(v2));

        PolygonA.Vertices[0] = Subject->Vertices[0];
        PolygonA.Vertices[1] = Subject->Vertices[1];
        PolygonA.Vertices[2] = Subject->Vertices[2];
        PolygonA.VertexCount = 3;

        // NOTE(babykaban): Constracting the overlap polygon
        polygon2 Poly = {};
        Poly.Vertices = (v2 *)malloc(SUBTRACTION_MAX_POINTS_PER_POLYGON*sizeof(v2));
        SutherlandHodgman(&PolygonA, Subtractor, &Poly);
        CleanUpOverlapPolygon(&Poly);

        f32 OverlappingArea = AbsoluteValue(PolygonSignedArea(&Poly));
        if((OverlappingArea > MinimalOverlapArea) && ((AbsoluteValue(AreaA) - OverlappingArea) > MinimalResultArea))
        {
            // NOTE(babykaban): Checks if any points of a subtractor lies on first vertex of subject if so rotate points in array
            HandleVertexCase(&PolygonA, &Poly);

            // NOTE(babykaban): Adding first vertices to the end, for algorithm to work properly
            PolygonA.Vertices[3] = PolygonA.Vertices[0];
            Poly.Vertices[Poly.VertexCount] = Poly.Vertices[0];
            ++PolygonA.VertexCount;
            ++Poly.VertexCount;

            // NOTE(babykaban): Find outside points for Subject polygon, because Overlap polygon will always be inside (on edge points are inside)
            s32 OutsidePointCount = 0;
            v2 *OutsidePoints = (v2 *)malloc(SUBTRACTION_MAX_POINTS_PER_POLYGON*sizeof(v2));
            FindOutsidePointsForF32(&PolygonA, &Poly, OutsidePoints, &OutsidePointCount);

            // NOTE(babykaban): Find intersection points for subject polygon
            s32 NewPointCount = 0;
            v2 *PointsToAdd = (v2 *)malloc(SUBTRACTION_MAX_POINTS_PER_POLYGON*sizeof(v2));
            points_between *PointTable = GetIntersectionPointsF32(&PolygonA, &Poly, PointsToAdd, &NewPointCount);

            // NOTE(babykaban): If no intersections was found record subtractor as a hole
            if(NewPointCount != 0)
            {
                // NOTE(babykaban): Insert points avoiding dublicats
                s32 PointIndex = 0;
                for(s32 I = 0;
                    I < 3;
                    ++I)
                {
                    points_between Between = PointTable[I];
                    while(Between.Count)
                    {
                        v2 p = PointsToAdd[PointIndex];
                        if(IsNewPoint(&PolygonA, p))
                        {
                            InsertPointBetween(&PolygonA, p, {Between.I0, Between.I1});

                            if(I == 0)
                            {
                                ++Between.I0;
                                ++Between.I1;
                                ++PointTable[1].I0;
                                ++PointTable[1].I1;
                                ++PointTable[2].I0;
                                ++PointTable[2].I1;
                            }
                            else if(I == 1)
                            {
                                ++Between.I0;
                                ++Between.I1;
                                ++PointTable[I + 1].I0;
                                ++PointTable[I + 1].I1;
                            }
                            else
                            {
                                ++Between.I0;
                                ++Between.I1;
                            }
                        }

                        Between.Count -= 1;
                        PointIndex += 1;
                    }
                }
                
                // NOTE(babykaban): Constract vertex table for both polygons
                vertex_info *PolygonAVertexInfo = (vertex_info *)malloc(PolygonA.VertexCount*sizeof(vertex_info));
                vertex_info *PolygonBVertexInfo = (vertex_info *)malloc(Poly.VertexCount*sizeof(vertex_info));

                BuildVertexTableForF32(&PolygonA, &Poly, OutsidePoints, OutsidePointCount, PolygonAVertexInfo);
                BuildVertexTableForF32(&Poly, &PolygonA, 0, 0, PolygonBVertexInfo);

                RecordCrossingPointsF32(PolygonAVertexInfo, PolygonA.VertexCount, PolygonBVertexInfo, Poly.VertexCount);

                // NOTE(babykaban): Construct resulting polygons
                Result.Success = ConstructPolygonsFor(&PolygonA, &Poly, PolygonAVertexInfo, PolygonBVertexInfo, &Result.Set);

                for(s32 PolygonIndex = 0;
                    PolygonIndex < Result.Set.PolygonCount;
                    ++PolygonIndex)
                {
                    polygon2 *P = Result.Set.Polygons + PolygonIndex;
                    RemoveDublicatPoints(P);

                    IsPolygonCollinearF32(P, CollinearValue);
                    if(P->VertexCount == 0)
                    {
                        for(s32 J = PolygonIndex;
                            J < (Result.Set.PolygonCount - 1);
                            ++J)
                        {
                            Result.Set.Polygons[J].VertexCount = Result.Set.Polygons[J + 1].VertexCount;
                            Copy(sizeof(v2)*SUBTRACTION_MAX_POINTS_PER_POLYGON, Result.Set.Polygons[J + 1].Vertices, Result.Set.Polygons[J].Vertices);
                        }

                        polygon2 *Remove = Result.Set.Polygons + (Result.Set.PolygonCount - 1);
                        free(Remove->Vertices);
                        --Result.Set.PolygonCount;
                        --PolygonIndex;
                    }
                }

                free(PolygonA.Vertices);
                free(Poly.Vertices);
                free(OutsidePoints);
                free(PointsToAdd);
                free(PolygonAVertexInfo);
                free(PolygonBVertexInfo);
            }
            else
            {
                free(PolygonA.Vertices);
                free(Poly.Vertices);
                free(OutsidePoints);
                free(PointsToAdd);

                Result.Set.Polygons[0].Vertices = (v2 *)malloc(sizeof(v2)*3);

                Result.Set.Polygons[0].Vertices[0] = Subject->Vertices[0];
                Result.Set.Polygons[0].Vertices[1] = Subject->Vertices[1];
                Result.Set.Polygons[0].Vertices[2] = Subject->Vertices[2];
                Result.Set.Polygons[0].VertexCount = 3;

                Result.Set.Polygons[0].HasHole = true;
                Result.Set.Polygons[0].HoleVertices = (v2 *)malloc(sizeof(v2)*3);
                Result.Set.Polygons[0].HoleVertices[0] = Subtractor->Vertices[0];
                Result.Set.Polygons[0].HoleVertices[1] = Subtractor->Vertices[1];
                Result.Set.Polygons[0].HoleVertices[2] = Subtractor->Vertices[2];
                Result.Set.Polygons[0].HoleVertexCount = 3;
                Result.Success = true;

                Result.Set.PolygonCount = 1;
            }
        }
        else
        {
            free(PolygonA.Vertices);
            free(Poly.Vertices);

            if((AbsoluteValue(AreaA) - OverlappingArea) < MinimalResultArea)
            {
                free(Result.Set.Polygons);
            }
            else
            {
                Result.Set.Polygons[0].Vertices = (v2 *)malloc(sizeof(v2)*3);
                Result.Set.Polygons[0].Vertices[0] = Subject->Vertices[0];
                Result.Set.Polygons[0].Vertices[1] = Subject->Vertices[1];
                Result.Set.Polygons[0].Vertices[2] = Subject->Vertices[2];
                Result.Set.Polygons[0].VertexCount = 3;
                Result.Set.PolygonCount = 1;
            }
            
            Result.Success = true;
        }
    }
    else
    {
        Result.Set.Polygons[0].Vertices = (v2 *)malloc(sizeof(v2)*3);
        Result.Set.Polygons[0].Vertices[0] = Subject->Vertices[0];
        Result.Set.Polygons[0].Vertices[1] = Subject->Vertices[1];
        Result.Set.Polygons[0].Vertices[2] = Subject->Vertices[2];
        Result.Set.Polygons[0].VertexCount = 3;

        Result.Set.PolygonCount = 1;
        Result.Success = true;
    }

    
    return(Result);
}

int main()
{
    triangle Subject = {};
    Subject.Vertices[0] = V2(-3.0f, 0.0f);
    Subject.Vertices[1] = V2(4.0f, 3.0f);
    Subject.Vertices[2] = V2(3.0f, -4.0f);

    triangle Subtractor = {};
    Subtractor.Vertices[0] = V2(-3.0f, -2.0f);
    Subtractor.Vertices[1] = V2(2.5f, 4.0f);
    Subtractor.Vertices[2] = V2(5.0f, -3.0f);

    subtract_result Result = SubtractTriangels(&Subject, &Subtractor, 6.25e-3f, 0.06f, 0.125f);
    if(Result.Success)
    {
        fprintf(stdout, "Result: \n");
        for(s32 PolygonIndex = 0;
            PolygonIndex < Result.Set.PolygonCount;
            ++PolygonIndex)
        {
            polygon2 *Poly = Result.Set.Polygons + PolygonIndex;
            fprintf(stdout, "%d. ", PolygonIndex);
            for(s32 VertexIndex = 0;
                VertexIndex < Poly->VertexCount;
                ++VertexIndex)
            {
                fprintf(stdout, "V2(%.2f, %.2f), ", Poly->Vertices[VertexIndex].x, Poly->Vertices[VertexIndex].y);
            }
            fprintf(stdout, "\n");
        }
    }
    else
    {
        fprintf(stdout, "fail");
    }
    
    return(0);
}
