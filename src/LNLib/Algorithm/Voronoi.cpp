/*
MapManager library for the conversion, manipulation and analysis
of maps used in Mobile Robotics research.
Copyright (C) 2005 Shane O'Sullivan

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

email: shaneosullivan1@gmail.com
*/

/*
 * Author:
 * 2024/04/06 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 
 *
 * Use of this source code is governed by a GPL-3.0 license that can be found in
 * the LICENSE file.
 */

#include "Voronoi.h"

#include <cmath>

LNLib::VoronoiDiagramGenerator::VoronoiDiagramGenerator()
{
	siteidx = 0;
	sites = 0;

	allMemoryList = new LNLib::FreeNodeArrayList;
	allMemoryList->memory = 0;
	allMemoryList->next = 0;
	currentMemoryBlock = allMemoryList;
	allEdges = 0;
	iteratorEdges = 0;
	sites = 0;
	ELhash = 0;
	PQhash = 0;

	minDistanceBetweenSites = 0;

	vertexLinks = 0;
	vertices = 0;
	finalVertexLinks = 0;
	finalVertices = 0;

	setGenerateVoronoi(true);
	setGenerateDelaunay(false);

	delaunayEdges = 0;
	iteratorDelaunayEdges = 0;
}

LNLib::VoronoiDiagramGenerator::~VoronoiDiagramGenerator()
{
	reset();
}

void LNLib::VoronoiDiagramGenerator::reset()
{
	cleanup();
	cleanupEdges();

	if (allMemoryList != 0)
		delete allMemoryList;

	if (finalVertices != 0)
		free(finalVertices);

	if (vertexLinks != 0)
		free(vertexLinks);

	if (vertices != 0)
		free(vertices);

	if (finalVertexLinks != 0)
		free(finalVertexLinks);

	allMemoryList = 0;
	finalVertices = 0;
	vertexLinks = 0;
	vertices = 0;
	finalVertexLinks = 0;
	triangles.clear();
}

std::vector<std::vector<int>> LNLib::VoronoiDiagramGenerator::getTriangles()
{
	return triangles;
}

void LNLib::VoronoiDiagramGenerator::setGenerateDelaunay(bool genDel)
{
	genDelaunay = genDel;
}

void LNLib::VoronoiDiagramGenerator::setGenerateVoronoi(bool genVor)
{
	genVoronoi = genVor;
}

bool LNLib::VoronoiDiagramGenerator::generateVoronoi(float* xValues, float* yValues, int numPoints, float minX,
	float maxX, float minY, float maxY, float minDist, bool genVertexInfo)
{
	cleanup();
	cleanupEdges();
	int i;

	minDistanceBetweenSites = minDist;

	nsites = numPoints;
	plot = 0;
	triangulate = 0;
	debug = 1;
	sorted = 0;
	if (sites != 0)
		free(sites);

	freeinit(&sfl, sizeof(LNLib::Site));

	sites = (struct Site*)myalloc(nsites * sizeof(*sites));

	if (finalVertices != 0)
		free(finalVertices);

	if (vertexLinks != 0)
		free(vertexLinks);

	if (vertices != 0)
		free(vertices);

	if (finalVertexLinks != 0)
		free(finalVertexLinks);

	vertexLinks = 0;
	vertices = 0;
	finalVertexLinks = 0;
	finalVertices = 0;

	sizeOfVertices = 0;
	sizeOfVertexLinks = 0;
	sizeOfFinalVertices = 0;
	sizeOfFinalVertexLinks = 0;

	if (sites == 0)
	{
		return false;
	}

	xmin = xValues[0];
	ymin = yValues[0];
	xmax = xValues[0];
	ymax = yValues[0];

	for (i = 0; i < nsites; i++)
	{
		sites[i].coord.x = xValues[i];
		sites[i].coord.y = yValues[i];
		sites[i].sitenbr = i;
		sites[i].refcnt = 0;

		if (xValues[i] < xmin)
			xmin = xValues[i];
		else if (xValues[i] > xmax)
			xmax = xValues[i];

		if (yValues[i] < ymin)
			ymin = yValues[i];
		else if (yValues[i] > ymax)
			ymax = yValues[i];

	}

	qsort(sites, nsites, sizeof(*sites), scomp); //undo

	siteidx = 0;
	geominit();
	float temp = 0;
	if (minX > maxX)
	{
		temp = minX;
		minX = maxX;
		maxX = temp;
	}
	if (minY > maxY)
	{
		temp = minY;
		minY = maxY;
		maxY = temp;
	}
	borderMinX = minX;
	borderMinY = minY;
	borderMaxX = maxX;
	borderMaxY = maxY;

	siteidx = 0;
	voronoi(genVertexInfo); //uncomment

	return true;
}

bool LNLib::VoronoiDiagramGenerator::ELinitialize()
{
	int i;
	freeinit(&hfl, sizeof * *ELhash);
	ELhashsize = 2 * sqrt_nsites;
	ELhash = (struct Halfedge**)myalloc(sizeof * ELhash * ELhashsize);

	if (ELhash == 0)
		return false;

	for (i = 0; i < ELhashsize; i += 1) ELhash[i] = (struct Halfedge*)NULL;
	ELleftend = HEcreate((struct Edge*)NULL, 0);
	ELrightend = HEcreate((struct Edge*)NULL, 0);
	ELleftend->ELleft = (struct Halfedge*)NULL;
	ELleftend->ELright = ELrightend;
	ELrightend->ELleft = ELleftend;
	ELrightend->ELright = (struct Halfedge*)NULL;
	ELhash[0] = ELleftend;
	ELhash[ELhashsize - 1] = ELrightend;

	return true;
}


struct LNLib::Halfedge* LNLib::VoronoiDiagramGenerator::HEcreate(struct Edge* e, int pm)
{
	struct Halfedge* answer;
	answer = (struct Halfedge*)getfree(&hfl);
	answer->ELedge = e;
	answer->ELpm = pm;
	answer->PQnext = (struct Halfedge*)NULL;
	answer->vertex = (struct Site*)NULL;
	answer->ELrefcnt = 0;
	return(answer);
}

void LNLib::VoronoiDiagramGenerator::ELinsert(struct	Halfedge* lb, struct Halfedge* newHe)
{
	newHe->ELleft = lb;
	newHe->ELright = lb->ELright;
	(lb->ELright)->ELleft = newHe;
	lb->ELright = newHe;
}

/* Get entry from hash table, pruning any deleted nodes */
struct LNLib::Halfedge* LNLib::VoronoiDiagramGenerator::ELgethash(int b)
{
	struct Halfedge* he;

	if (b < 0 || b >= ELhashsize)
		return((struct Halfedge*)NULL);
	he = ELhash[b];
	if (he == (struct Halfedge*)NULL || he->ELedge != (struct Edge*)DELETED)
		return (he);

	/* Hash table points to deleted half edge.  Patch as necessary. */
	ELhash[b] = (struct Halfedge*)NULL;
	if ((he->ELrefcnt -= 1) == 0)
		makefree((Freenode*)he, &hfl);
	return ((struct Halfedge*)NULL);
}

struct LNLib::Halfedge* LNLib::VoronoiDiagramGenerator::ELleftbnd(struct PointVDG* p)
{
	int i, bucket;
	struct Halfedge* he;

	/* Use hash table to get close to desired halfedge */
	bucket = (int)((p->x - xmin) / deltax * ELhashsize);	//use the hash function to find the place in the hash map that this HalfEdge should be

	if (bucket < 0) bucket = 0;					//make sure that the bucket position in within the range of the hash array
	if (bucket >= ELhashsize) bucket = ELhashsize - 1;

	he = ELgethash(bucket);
	if (he == (struct Halfedge*)NULL)			//if the HE isn't found, search backwards and forwards in the hash map for the first non-null entry
	{
		for (i = 1; 1; i += 1)
		{
			if ((he = ELgethash(bucket - i)) != (struct Halfedge*)NULL)
				break;
			if ((he = ELgethash(bucket + i)) != (struct Halfedge*)NULL)
				break;
		};
		totalsearch += i;
	};
	ntry += 1;
	/* Now search linear list of halfedges for the correct one */
	if (he == ELleftend || (he != ELrightend && right_of(he, p)))
	{
		do
		{
			he = he->ELright;
		} while (he != ELrightend && right_of(he, p));	//keep going right on the list until either the end is reached, or you find the 1st edge which the point
		he = he->ELleft;				//isn't to the right of
	}
	else 							//if the point is to the left of the HalfEdge, then search left for the HE just to the left of the point
		do
		{
			he = he->ELleft;
		} while (he != ELleftend && !right_of(he, p));

		/* Update hash table and reference counts */
		if (bucket > 0 && bucket < ELhashsize - 1)
		{
			if (ELhash[bucket] != (struct Halfedge*)NULL)
			{
				ELhash[bucket]->ELrefcnt -= 1;
			}
			ELhash[bucket] = he;
			ELhash[bucket]->ELrefcnt += 1;
		};
		return (he);
}


/* This delete routine can't reclaim node, since pointers from hash
table may be present.   */
void LNLib::VoronoiDiagramGenerator::ELdelete(struct Halfedge* he)
{
	(he->ELleft)->ELright = he->ELright;
	(he->ELright)->ELleft = he->ELleft;
	he->ELedge = (struct Edge*)DELETED;
}


struct LNLib::Halfedge* LNLib::VoronoiDiagramGenerator::ELright(struct Halfedge* he)
{
	return (he->ELright);
}

struct LNLib::Halfedge* LNLib::VoronoiDiagramGenerator::ELleft(struct Halfedge* he)
{
	return (he->ELleft);
}

struct LNLib::Site* LNLib::VoronoiDiagramGenerator::leftreg(struct Halfedge* he)
{
	if (he->ELedge == (struct Edge*)NULL)
		return(bottomsite);
	return(he->ELpm == le ?
		he->ELedge->reg[le] : he->ELedge->reg[re]);
}

struct LNLib::Site* LNLib::VoronoiDiagramGenerator::rightreg(struct Halfedge* he)
{
	if (he->ELedge == (struct Edge*)NULL) //if this halfedge has no edge, return the bottom site (whatever that is)
		return(bottomsite);

	//if the ELpm field is zero, return the site 0 that this edge bisects, otherwise return site number 1
	return(he->ELpm == le ? he->ELedge->reg[re] : he->ELedge->reg[le]);
}

void LNLib::VoronoiDiagramGenerator::geominit()
{
	float sn;

	freeinit(&efl, sizeof(Edge));
	nvertices = 0;
	nedges = 0;
	sn = (float)nsites + 4;
	sqrt_nsites = (int)sqrt(sn);
	deltay = ymax - ymin;
	deltax = xmax - xmin;
}

struct LNLib::Edge* LNLib::VoronoiDiagramGenerator::bisect(struct Site* s1, struct Site* s2)
{
	float dx, dy, adx, ady;
	struct Edge* newedge;

	newedge = (struct Edge*)getfree(&efl);

	newedge->reg[0] = s1; //store the sites that this edge is bisecting
	newedge->reg[1] = s2;
	ref(s1);
	ref(s2);
	newedge->ep[0] = (struct Site*)NULL; //to begin with, there are no endpoints on the bisector - it goes to infinity
	newedge->ep[1] = (struct Site*)NULL;

	dx = s2->coord.x - s1->coord.x;			//get the difference in x dist between the sites
	dy = s2->coord.y - s1->coord.y;
	adx = dx > 0 ? dx : -dx;					//make sure that the difference in positive
	ady = dy > 0 ? dy : -dy;
	newedge->c = (float)(s1->coord.x * dx + s1->coord.y * dy + (dx * dx + dy * dy) * 0.5);//get the slope of the line

	if (adx > ady)
	{
		newedge->a = 1.0; newedge->b = dy / dx; newedge->c /= dx;//set formula of line, with x fixed to 1
	}
	else
	{
		newedge->b = 1.0; newedge->a = dx / dy; newedge->c /= dy;//set formula of line, with y fixed to 1
	};

	newedge->edgenbr = nedges;
	out_bisector(newedge);
	nedges += 1;
	return(newedge);
}

//create a new site where the HalfEdges el1 and el2 intersect - note that the PointVDG in the argument list is not used, don't know why it's there
struct LNLib::Site* LNLib::VoronoiDiagramGenerator::intersect(struct Halfedge* el1, struct Halfedge* el2, struct PointVDG* p)
{
	struct	Edge* e1, * e2, * e;
	struct  Halfedge* el;
	float d, xint, yint;
	int right_of_site;
	struct Site* v;

	e1 = el1->ELedge;
	e2 = el2->ELedge;
	if (e1 == (struct Edge*)NULL || e2 == (struct Edge*)NULL)
		return ((struct Site*)NULL);

	//if the two edges bisect the same parent, return null
	if (e1->reg[1] == e2->reg[1])
		return ((struct Site*)NULL);

	d = e1->a * e2->b - e1->b * e2->a;
	if (-1.0e-10 < d && d < 1.0e-10)
		return ((struct Site*)NULL);

	xint = (e1->c * e2->b - e2->c * e1->b) / d;
	yint = (e2->c * e1->a - e1->c * e2->a) / d;

	if ((e1->reg[1]->coord.y < e2->reg[1]->coord.y) ||
		(e1->reg[1]->coord.y == e2->reg[1]->coord.y &&
			e1->reg[1]->coord.x < e2->reg[1]->coord.x))
	{
		el = el1;
		e = e1;
	}
	else
	{
		el = el2;
		e = e2;
	};

	right_of_site = xint >= e->reg[1]->coord.x;
	if ((right_of_site && el->ELpm == le) || (!right_of_site && el->ELpm == re))
		return ((struct Site*)NULL);

	//create a new site at the point of intersection - this is a new vector event waiting to happen
	v = (struct Site*)getfree(&sfl);
	v->refcnt = 0;
	v->coord.x = xint;
	v->coord.y = yint;
	return(v);
}

/* returns 1 if p is to right of halfedge e */
int LNLib::VoronoiDiagramGenerator::right_of(struct Halfedge* el, struct PointVDG* p)
{
	struct Edge* e;
	struct Site* topsite;
	int right_of_site, above, fast;
	float dxp, dyp, dxs, t1, t2, t3, yl;

	e = el->ELedge;
	topsite = e->reg[1];
	right_of_site = p->x > topsite->coord.x;
	if (right_of_site && el->ELpm == le) return(1);
	if (!right_of_site && el->ELpm == re) return (0);

	if (e->a == 1.0)
	{
		dyp = p->y - topsite->coord.y;
		dxp = p->x - topsite->coord.x;
		fast = 0;
		if ((!right_of_site & (e->b < 0.0)) | (right_of_site & (e->b >= 0.0)))
		{
			above = dyp >= e->b * dxp;
			fast = above;
		}
		else
		{
			above = p->x + p->y * e->b > e->c;
			if (e->b < 0.0) above = !above;
			if (!above) fast = 1;
		};
		if (!fast)
		{
			dxs = topsite->coord.x - (e->reg[0])->coord.x;
			above = e->b * (dxp * dxp - dyp * dyp) <
				dxs * dyp * (1.0 + 2.0 * dxp / dxs + e->b * e->b);
			if (e->b < 0.0) above = !above;
		};
	}
	else  /*e->b==1.0 */
	{
		yl = e->c - e->a * p->x;
		t1 = p->y - yl;
		t2 = p->x - topsite->coord.x;
		t3 = yl - topsite->coord.y;
		above = t1 * t1 > t2 * t2 + t3 * t3;
	};
	return (el->ELpm == le ? above : !above);
}

void LNLib::VoronoiDiagramGenerator::endpoint(struct Edge* e, int lr, struct Site* s)
{
	e->ep[lr] = s;
	ref(s);
	if (e->ep[re - lr] == (struct Site*)NULL)
		return;

	clip_line(e);

	deref(e->reg[le]);
	deref(e->reg[re]);
	makefree((Freenode*)e, &efl);
}

float LNLib::VoronoiDiagramGenerator::dist(struct Site* s, struct Site* t)
{
	float dx, dy;
	dx = s->coord.x - t->coord.x;
	dy = s->coord.y - t->coord.y;
	return (float)(sqrt(dx * dx + dy * dy));
}

void LNLib::VoronoiDiagramGenerator::makevertex(struct Site* v)
{
	v->sitenbr = nvertices;
	insertVertexAddress(nvertices, v);
	nvertices += 1;
}

void LNLib::VoronoiDiagramGenerator::out_triple(Site* s1, Site* s2, Site* s3)
{
	std::vector<int> v(3);
	v[0] = s1->sitenbr;
	v[1] = s2->sitenbr;
	v[2] = s3->sitenbr;
	triangles.emplace_back(v);
}

void LNLib::VoronoiDiagramGenerator::deref(struct Site* v)
{
	v->refcnt -= 1;
	if (v->refcnt == 0)
		makefree((Freenode*)v, &sfl);
}

void LNLib::VoronoiDiagramGenerator::ref(struct Site* v)
{
	v->refcnt += 1;
	v->overallRefcnt += 1;
}

//push the HalfEdge into the ordered linked list of vertices
void LNLib::VoronoiDiagramGenerator::PQinsert(struct Halfedge* he, struct Site* v, float offset)
{
	struct Halfedge* last, * next;

	he->vertex = v;
	ref(v);
	he->ystar = (float)(v->coord.y + offset);
	last = &PQhash[PQbucket(he)];
	while ((next = last->PQnext) != (struct Halfedge*)NULL &&
		(he->ystar > next->ystar ||
			(he->ystar == next->ystar && v->coord.x > next->vertex->coord.x)))
	{
		last = next;
	};
	he->PQnext = last->PQnext;
	last->PQnext = he;
	PQcount += 1;
}

//remove the HalfEdge from the list of vertices 
void LNLib::VoronoiDiagramGenerator::PQdelete(struct Halfedge* he)
{
	struct Halfedge* last;

	if (he->vertex != (struct Site*)NULL)
	{
		last = &PQhash[PQbucket(he)];
		while (last->PQnext != he)
			last = last->PQnext;

		last->PQnext = he->PQnext;
		PQcount -= 1;
		deref(he->vertex);
		he->vertex = (struct Site*)NULL;
	};
}

int LNLib::VoronoiDiagramGenerator::PQbucket(struct Halfedge* he)
{
	int bucket;

	bucket = (int)((he->ystar - ymin) / deltay * PQhashsize);
	if (bucket < 0) bucket = 0;
	if (bucket >= PQhashsize) bucket = PQhashsize - 1;
	if (bucket < PQmin) PQmin = bucket;
	return(bucket);
}

int LNLib::VoronoiDiagramGenerator::PQempty()
{
	return(PQcount == 0);
}

struct LNLib::PointVDG LNLib::VoronoiDiagramGenerator::PQ_min()
{
	struct PointVDG answer;

	while (PQhash[PQmin].PQnext == (struct Halfedge*)NULL) { PQmin += 1; };
	answer.x = PQhash[PQmin].PQnext->vertex->coord.x;
	answer.y = PQhash[PQmin].PQnext->ystar;
	return (answer);
}

struct LNLib::Halfedge* LNLib::VoronoiDiagramGenerator::PQextractmin()
{
	struct Halfedge* curr;

	curr = PQhash[PQmin].PQnext;
	PQhash[PQmin].PQnext = curr->PQnext;
	PQcount -= 1;
	return(curr);
}

bool LNLib::VoronoiDiagramGenerator::PQinitialize()
{
	int i;

	PQcount = 0;
	PQmin = 0;
	PQhashsize = 4 * sqrt_nsites;
	PQhash = (struct Halfedge*)myalloc(PQhashsize * sizeof * PQhash);

	if (PQhash == 0)
		return false;

	for (i = 0; i < PQhashsize; i += 1) PQhash[i].PQnext = (struct Halfedge*)NULL;

	return true;
}

void LNLib::VoronoiDiagramGenerator::freeinit(struct Freelist* fl, int size)
{
	fl->head = (struct Freenode*)NULL;
	fl->nodesize = size;
}

char* LNLib::VoronoiDiagramGenerator::getfree(struct Freelist* fl)
{
	int i;
	struct Freenode* t;

	if (fl->head == (struct Freenode*)NULL)
	{
		t = (struct Freenode*)myalloc(sqrt_nsites * fl->nodesize);

		if (t == 0)
			return 0;

		currentMemoryBlock->next = new FreeNodeArrayList;
		currentMemoryBlock = currentMemoryBlock->next;
		currentMemoryBlock->memory = t;
		currentMemoryBlock->next = 0;

		for (i = 0; i < sqrt_nsites; i += 1)
			makefree((struct Freenode*)((char*)t + i * fl->nodesize), fl);
	};
	t = fl->head;
	fl->head = (fl->head)->nextfree;
	return((char*)t);
}

void LNLib::VoronoiDiagramGenerator::makefree(struct Freenode* curr, struct Freelist* fl)
{
	curr->nextfree = fl->head;
	fl->head = curr;
}

void LNLib::VoronoiDiagramGenerator::cleanup()
{
	if (sites != 0)
	{
		free(sites);
		sites = 0;
	}

	FreeNodeArrayList* current = 0, * prev = 0;

	current = prev = allMemoryList;

	if (current != 0)
	{
		while (current->next != 0)
		{
			prev = current;
			current = current->next;
			free(prev->memory);
			delete prev;
			prev = 0;
		}
	}
	allMemoryList = 0;

	if (current != 0 && current->memory != 0)
	{
		free(current->memory);
		delete current;
	}

	allMemoryList = new FreeNodeArrayList;
	allMemoryList->next = 0;
	allMemoryList->memory = 0;
	currentMemoryBlock = allMemoryList;

	if (ELhash != 0)
	{
		free(ELhash);
		ELhash = 0;
	}

	if (PQhash != 0)
	{
		free(PQhash);
		PQhash = 0;
	}
}

void LNLib::VoronoiDiagramGenerator::cleanupEdges()
{
	GraphEdge* geCurrent = 0, * gePrev = 0;
	geCurrent = gePrev = allEdges;

	while (geCurrent != 0 && geCurrent->next != 0)
	{
		gePrev = geCurrent;
		geCurrent = geCurrent->next;
		delete gePrev;
	}

	allEdges = 0;

	geCurrent = gePrev = delaunayEdges;

	while (geCurrent != 0 && geCurrent->next != 0)
	{
		gePrev = geCurrent;
		geCurrent = geCurrent->next;
		delete gePrev;
	}

	delaunayEdges = 0;
}

void LNLib::VoronoiDiagramGenerator::pushGraphEdge(float x1, float y1, float x2, float y2)
{
	if (genVoronoi)
	{
		GraphEdge* newEdge = new GraphEdge;
		newEdge->next = allEdges;
		allEdges = newEdge;
		newEdge->x1 = x1;
		newEdge->y1 = y1;
		newEdge->x2 = x2;
		newEdge->y2 = y2;
	}
}

void LNLib::VoronoiDiagramGenerator::pushDelaunayGraphEdge(float x1, float y1, float x2, float y2)
{
	if (std::sqrt(((x2 - x1) * (x2 - x1)) + ((y2 - y1) * (y2 - y1))) < minDistanceBetweenSites)
	{
		return;
	}
	GraphEdge* newEdge = new GraphEdge;
	newEdge->next = delaunayEdges;
	delaunayEdges = newEdge;
	newEdge->x1 = x1;
	newEdge->y1 = y1;
	newEdge->x2 = x2;
	newEdge->y2 = y2;
}

char* LNLib::VoronoiDiagramGenerator::myalloc(unsigned n)
{
	char* t = 0;
	t = (char*)malloc(n);
	total_alloc += n;
	return(t);
}

/* for those who don't have Cherry's plot */
/* #include <plot.h> */
void LNLib::VoronoiDiagramGenerator::line(float x1, float y1, float x2, float y2)
{
	pushGraphEdge(x1, y1, x2, y2);

}
void LNLib::VoronoiDiagramGenerator::circle(float x, float y, float radius) {}

void LNLib::VoronoiDiagramGenerator::range(float minX, float minY, float maxX, float maxY) {}

void LNLib::VoronoiDiagramGenerator::out_bisector(struct Edge* e)
{
	if (genDelaunay)
	{
		pushDelaunayGraphEdge(e->reg[0]->coord.x, e->reg[0]->coord.y,
			e->reg[1]->coord.x, e->reg[1]->coord.y);
	}
}

void LNLib::VoronoiDiagramGenerator::plotinit()
{
	float dx, dy, d;

	dy = ymax - ymin;
	dx = xmax - xmin;
	d = (float)((dx > dy ? dx : dy) * 1.1);
	pxmin = (float)(xmin - (d - dx) / 2.0);
	pxmax = (float)(xmax + (d - dx) / 2.0);
	pymin = (float)(ymin - (d - dy) / 2.0);
	pymax = (float)(ymax + (d - dy) / 2.0);
	cradius = (float)((pxmax - pxmin) / 350.0);
	range(pxmin, pymin, pxmax, pymax);
}

void LNLib::VoronoiDiagramGenerator::clip_line(struct Edge* e)
{
	struct Site* s1, * s2;
	float x1 = 0, x2 = 0, y1 = 0, y2 = 0, temp = 0;
	Site* v1 = 0, * v2 = 0;
	bool needNewVertex1 = false, needNewVertex2 = false;

	x1 = e->reg[0]->coord.x;
	x2 = e->reg[1]->coord.x;
	y1 = e->reg[0]->coord.y;
	y2 = e->reg[1]->coord.y;

	//if the distance between the two points this line was created from is less than 
	//the square root of 2, then ignore it
	if (std::sqrt(((x2 - x1) * (x2 - x1)) + ((y2 - y1) * (y2 - y1))) < minDistanceBetweenSites)
	{
		return;
	}
	pxmin = borderMinX;
	pxmax = borderMaxX;
	pymin = borderMinY;
	pymax = borderMaxY;

	//	printf("\nEdge (%d), minX = %f, maxX = %f, minY = %f, maxY = %f",e->edgenbr,pxmin, pxmax, pymin,pymax);
	if (e->a == 1.0 && e->b >= 0.0)
	{
		s1 = e->ep[1];
		s2 = e->ep[0];
	}
	else
	{
		s1 = e->ep[0];
		s2 = e->ep[1];
	};

	v1 = s1;
	v2 = s2;

	if (e->a == 1.0)
	{
		y1 = pymin;
		if (s1 != (struct Site*)NULL && s1->coord.y > pymin)
		{
			y1 = s1->coord.y;
		}
		else
		{
			needNewVertex1 = true;
		}

		if (y1 > pymax)
		{
			y1 = pymax;
			needNewVertex1 = true;
		}

		x1 = e->c - e->b * y1;
		y2 = pymax;
		if (s2 != (struct Site*)NULL && s2->coord.y < pymax)
		{
			y2 = s2->coord.y;
		}
		else
		{
			needNewVertex2 = true;
		}


		if (y2 < pymin)
		{
			y2 = pymin;
			needNewVertex2 = true;
		}

		x2 = (e->c) - (e->b) * y2;
		if (((x1 > pxmax) & (x2 > pxmax)) | ((x1 < pxmin) & (x2 < pxmin)))
		{
			return;
		}
		if (x1 > pxmax)
		{
			x1 = pxmax; y1 = (e->c - x1) / e->b; needNewVertex1 = true;
		}
		if (x1 < pxmin)
		{
			x1 = pxmin; y1 = (e->c - x1) / e->b; needNewVertex1 = true;
		}
		if (x2 > pxmax)
		{
			x2 = pxmax; y2 = (e->c - x2) / e->b; needNewVertex2 = true;
		}
		if (x2 < pxmin)
		{
			x2 = pxmin; y2 = (e->c - x2) / e->b; needNewVertex2 = true;
		}
	}
	else
	{
		x1 = pxmin;
		if (s1 != (struct Site*)NULL && s1->coord.x > pxmin)
		{
			x1 = s1->coord.x;
		}
		else
		{
			needNewVertex1 = true;
		}

		if (x1 > pxmax)
		{
			x1 = pxmax;
			needNewVertex1 = true;
		}
		y1 = e->c - e->a * x1;
		x2 = pxmax;
		if (s2 != (struct Site*)NULL && s2->coord.x < pxmax)
		{
			x2 = s2->coord.x;
		}
		else
		{
			needNewVertex2 = true;
		}

		if (x2 < pxmin)
		{
			x2 = pxmin;
			needNewVertex2 = true;
		}
		y2 = e->c - e->a * x2;
		if (((y1 > pymax) & (y2 > pymax)) | ((y1 < pymin) & (y2 < pymin)))
		{
			return;
		}
		if (y1 > pymax)
		{
			y1 = pymax; x1 = (e->c - y1) / e->a; needNewVertex1 = true;
		}
		if (y1 < pymin)
		{
			y1 = pymin; x1 = (e->c - y1) / e->a; needNewVertex1 = true;
		}
		if (y2 > pymax)
		{
			y2 = pymax; x2 = (e->c - y2) / e->a; needNewVertex2 = true;
		}
		if (y2 < pymin)
		{
			y2 = pymin; x2 = (e->c - y2) / e->a; needNewVertex2 = true;
		}
	}

	if (!((x1 == x2 && x2 == pxmin) || (x1 == x2 && x2 == pxmax) ||
		(y1 == y2 && y2 == pymin) || (y1 == y2 && y2 == pymax)))
	{
		pushGraphEdge(x1, y1, x2, y2);
		if (needNewVertex1)
		{
			//printf("\nCreate new vertex 1 
			v1 = (struct Site*)getfree(&sfl);
			v1->refcnt = 0;
			v1->coord.x = x1;
			v1->coord.y = y1;
			makevertex(v1);

		}

		if (needNewVertex2)
		{
			v2 = (struct Site*)getfree(&sfl);
			v2->refcnt = 0;
			v2->coord.x = x2;
			v2->coord.y = y2;
			makevertex(v2);
		}
		insertVertexLink(v1->sitenbr, v2->sitenbr);
	}
}

/* implicit parameters: nsites, sqrt_nsites, xmin, xmax, ymin, ymax,
deltax, deltay (can all be estimates).
Performance suffers if they are wrong; better to make nsites,
deltax, and deltay too big than too small.  (?) */

bool LNLib::VoronoiDiagramGenerator::voronoi(bool genVertexInfo)
{
	struct Site* newsite, * bot, * top, * temp, * p;
	struct Site* v;
	struct PointVDG newintstar;
	int pm;
	struct Halfedge* lbnd, * rbnd, * llbnd, * rrbnd, * bisector;
	struct Edge* e;

	PQinitialize();
	bottomsite = nextone();
	bool retval = ELinitialize();

	if (!retval)
	{
		return false;
	}

	newsite = nextone();

	long counter = 0;
	while (1)
	{
		if (!PQempty())
			newintstar = PQ_min();

		//if the lowest site has a smaller y value than the lowest vector intersection, process the site
		//otherwise process the vector intersection		

		if (newsite != (struct Site*)NULL && (PQempty() || newsite->coord.y < newintstar.y
			|| (newsite->coord.y == newintstar.y && newsite->coord.x < newintstar.x)))
		{/* new site is smallest - this is a site event*/
		//	out_site(newsite);						//output the site
			lbnd = ELleftbnd(&(newsite->coord));				//get the first HalfEdge to the LEFT of the new site
			rbnd = ELright(lbnd);						//get the first HalfEdge to the RIGHT of the new site
			bot = rightreg(lbnd);						//if this halfedge has no edge, , bot = bottom site (whatever that is)
			e = bisect(bot, newsite);					//create a new edge that bisects 

			bisector = HEcreate(e, le);					//create a new HalfEdge, setting its ELpm field to 0			
			ELinsert(lbnd, bisector);					//insert this new bisector edge between the left and right vectors in a linked list	

			if ((p = intersect(lbnd, bisector)) != (struct Site*)NULL) 	//if the new bisector intersects with the left edge, remove the left edge's vertex, and put in the new one
			{
				PQdelete(lbnd);
				PQinsert(lbnd, p, dist(p, newsite));
			};
			lbnd = bisector;
			bisector = HEcreate(e, re);					//create a new HalfEdge, setting its ELpm field to 1
			ELinsert(lbnd, bisector);					//insert the new HE to the right of the original bisector earlier in the IF stmt

			if ((p = intersect(bisector, rbnd)) != (struct Site*)NULL)	//if this new bisector intersects with the
			{
				PQinsert(bisector, p, dist(p, newsite));			//push the HE into the ordered linked list of vertices
			};

			newsite = nextone();
		}
		else if (!PQempty()) /* intersection is smallest - this is a vector event */
		{
			lbnd = PQextractmin();						//pop the HalfEdge with the lowest vector off the ordered list of vectors				
			llbnd = ELleft(lbnd);						//get the HalfEdge to the left of the above HE
			rbnd = ELright(lbnd);						//get the HalfEdge to the right of the above HE
			rrbnd = ELright(rbnd);						//get the HalfEdge to the right of the HE to the right of the lowest HE 
			bot = leftreg(lbnd);						//get the Site to the left of the left HE which it bisects
			top = rightreg(rbnd);						//get the Site to the right of the right HE which it bisects

			out_triple(bot, top, rightreg(lbnd));

			v = lbnd->vertex;						//get the vertex that caused this event
			makevertex(v);							//set the vertex number - couldn't do this earlier since we didn't know when it would be processed
			endpoint(lbnd->ELedge, lbnd->ELpm, v);	//set the endpoint of the left HalfEdge to be this vector
			endpoint(rbnd->ELedge, rbnd->ELpm, v);	//set the endpoint of the right HalfEdge to be this vector
			ELdelete(lbnd);							//mark the lowest HE for deletion - can't delete yet because there might be pointers to it in Hash Map	
			PQdelete(rbnd);							//remove all vertex events to do with the  right HE
			ELdelete(rbnd);							//mark the right HE for deletion - can't delete yet because there might be pointers to it in Hash Map	
			pm = le;								//set the pm variable to zero

			if (bot->coord.y > top->coord.y)		//if the site to the left of the event is higher than the Site
			{										//to the right of it, then swap them and set the 'pm' variable to 1
				temp = bot;
				bot = top;
				top = temp;
				pm = re;
			}
			e = bisect(bot, top);					//create an Edge (or line) that is between the two Sites. This creates
			//the formula of the line, and assigns a line number to it
			bisector = HEcreate(e, pm);				//create a HE from the Edge 'e', and make it point to that edge with its ELedge field
			ELinsert(llbnd, bisector);				//insert the new bisector to the right of the left HE
			endpoint(e, re - pm, v);					//set one endpoint to the new edge to be the vector point 'v'.
			//If the site to the left of this bisector is higher than the right
			//Site, then this endpoint is put in position 0; otherwise in pos 1
			deref(v);								//delete the vector 'v'

			//if left HE and the new bisector don't intersect, then delete the left HE, and reinsert it 
			if ((p = intersect(llbnd, bisector)) != (struct Site*)NULL)
			{
				PQdelete(llbnd);
				PQinsert(llbnd, p, dist(p, bot));
			};

			//if right HE and the new bisector don't intersect, then reinsert it 
			if ((p = intersect(bisector, rrbnd)) != (struct Site*)NULL)
			{
				PQinsert(bisector, p, dist(p, bot));
			};
		}
		else break;
	};

	for (lbnd = ELright(ELleftend); lbnd != ELrightend; lbnd = ELright(lbnd))
	{
		e = lbnd->ELedge;

		clip_line(e);
	};

	//count the total number of 
	long i = 0;

	if (genVertexInfo)
	{
		counter = 0;
		sizeOfFinalVertices = 0;
		for (i = 0; i < sizeOfVertexLinks; i++)
		{
			if (vertexLinks[i].count == 1 || vertexLinks[i].count == 3)
			{
				sizeOfFinalVertices++;
			}
		}

		finalVertices = (PointVDG*)myalloc(sizeOfFinalVertices * sizeof(PointVDG));
		counter = 0;
		for (i = 0; i < sizeOfVertexLinks; i++)
		{
			if (vertexLinks[i].count == 1 || vertexLinks[i].count == 3 && vertices[i] != 0)
			{
				finalVertices[counter] = vertices[i]->coord;
				counter++;
			}
		}
		generateVertexLinks();

	}
	cleanup();
	return true;
}

void LNLib::VoronoiDiagramGenerator::insertVertexAddress(long vertexNum, struct Site* address)
{
	if (vertices == 0)
	{
		vertices = (struct Site**)myalloc(4000 * sizeof(*vertices));
		sizeOfVertices = 4000;
		for (int i = sizeOfVertices - 4000; i < sizeOfVertices; i++)
		{
			vertices[i] = 0;
		}
	}
	//if the site address being entered is past the end of the array, then grow the array
	while (vertexNum > sizeOfVertices - 1)
	{
		vertices = (Site**)realloc(vertices, (sizeOfVertices + 4000) * sizeof(Site*));
		sizeOfVertices += 4000;

		for (int i = sizeOfVertices - 4000; i < sizeOfVertices; i++)
		{
			vertices[i] = 0;
		}
	}
	vertices[vertexNum] = address;
}

void LNLib::VoronoiDiagramGenerator::insertVertexLink(long vertexNum, long vertexLinkedTo)
{
	if (vertexLinks == 0)
	{
		vertexLinks = (struct Point3*)myalloc(4000 * sizeof(*vertexLinks));
		sizeOfVertexLinks = 4000;
		for (int i = 0; i < sizeOfVertexLinks; i++)
		{
			vertexLinks[i].x = vertexLinks[i].y = vertexLinks[i].z = -1; //initialise all elements in the array to -1
			vertexLinks[i].count = 0;
		}
	}
	//if the site address being entered is past the end of the array, then grow the array
	while (vertexNum > sizeOfVertexLinks - 1 || vertexLinkedTo > sizeOfVertexLinks - 1)
	{
		vertexLinks = (Point3*)realloc(vertexLinks, (sizeOfVertexLinks + 4000) * sizeof(Point3));
		sizeOfVertexLinks += 4000;

		for (int i = sizeOfVertexLinks - 4000; i < sizeOfVertexLinks; i++)
		{
			vertexLinks[i].x = vertexLinks[i].y = vertexLinks[i].z = -1; //initialise all elements in the array to -1
			vertexLinks[i].count = 0;
		}
	}

	if (vertexLinks[vertexNum].x == -1)
	{
		vertexLinks[vertexNum].x = (float)vertexLinkedTo;
		vertexLinks[vertexNum].count = 1;
	}
	else if (vertexLinks[vertexNum].y == -1)
	{
		vertexLinks[vertexNum].y = (float)vertexLinkedTo;
		vertexLinks[vertexNum].count = 2;
	}
	else if (vertexLinks[vertexNum].z == -1)
	{
		vertexLinks[vertexNum].z = (float)vertexLinkedTo;
		vertexLinks[vertexNum].count = 3;
	}

	if (vertexLinks[vertexLinkedTo].x == -1)
	{
		vertexLinks[vertexLinkedTo].x = (float)vertexNum;
		vertexLinks[vertexLinkedTo].count = 1;
	}
	else if (vertexLinks[vertexLinkedTo].y == -1)
	{
		vertexLinks[vertexLinkedTo].y = (float)vertexNum;
		vertexLinks[vertexLinkedTo].count = 2;
	}
	else if (vertexLinks[vertexLinkedTo].z == -1)
	{
		vertexLinks[vertexLinkedTo].z = (float)vertexNum;
		vertexLinks[vertexLinkedTo].count = 3;
	}
}

void LNLib::VoronoiDiagramGenerator::generateVertexLinks()
{
	long i = 0, j = 0;
	if (finalVertexLinks != 0)
	{
		free(finalVertexLinks);
		finalVertexLinks = 0;
	}

	if (vertices == 0)
	{
		return;
	}

	sizeOfFinalVertexLinks = (sizeOfVertices > sizeOfVertexLinks) ? sizeOfVertices : sizeOfVertexLinks;
	finalVertexLinks = (struct LNLib::VertexLink*)myalloc(sizeOfFinalVertexLinks * sizeof(LNLib::VertexLink));

	if (finalVertexLinks == 0)
		return;
	for (i = 0; i < sizeOfFinalVertexLinks; i++)
	{
		finalVertexLinks[i].count = 0;
	}

	long* count1vertices = (long*)myalloc(sizeOfFinalVertexLinks * sizeof(long));
	long* count3vertices = (long*)myalloc(sizeOfFinalVertexLinks * sizeof(long));

	//if we couldn't get the memory we need, return
	if (count1vertices == 0 || count3vertices == 0)
	{
		if (count3vertices != 0)
		{
			delete[] count3vertices;
		}
		if (count1vertices != 0)
		{
			delete[] count1vertices;
		}
		return;
	}

	//initialise the two arrays
	for (i = 0; i < sizeOfFinalVertexLinks; i++)
	{
		count1vertices[i] = 0;
		count3vertices[i] = 0;
	}
	long count1counter = 0, count3counter = 0;
	long count3links[3] = { -1 };
	for (i = 0; i < sizeOfVertices; i++)
	{
		if (vertices[i] == 0)
			break;

		finalVertexLinks[i].coord = vertices[i]->coord;
	}

	for (i = 0; i < sizeOfVertexLinks; i++)
	{
		switch (vertexLinks[i].count)
		{
		case 1:count1vertices[count1counter] = i;
			count1counter++;
			break;
		case 3:count3vertices[count3counter] = i;
			count3counter++;
			break;
		}
	}
	//first, we go from all leaf nodes, those with just one edge, to either the next leaf or the
	//next node with 3 edges
	long currentVertex = 0, prevVertex = 0;

	for (i = 0; i < count1counter; i++)
	{
		if (vertexLinks[count1vertices[i]].count != 1)//if we've already been here, ignore this vertex
			continue;
		//	LOG<<i;

		currentVertex = (long)vertexLinks[count1vertices[i]].x;//use the x variable, since it should be the only one set
		prevVertex = count1vertices[i];
		vertexLinks[count1vertices[i]].count--;	//don't revisit this site

		while (currentVertex != -1 && currentVertex < sizeOfVertexLinks &&
			vertexLinks[currentVertex].count != 1 && vertexLinks[currentVertex].count != 3)
		{
			if (vertexLinks[currentVertex].count == 2)
			{
				if (vertexLinks[currentVertex].x == prevVertex)
				{
					vertexLinks[currentVertex].x = vertexLinks[currentVertex].y;
					vertexLinks[currentVertex].y = -1;
					vertexLinks[currentVertex].count--;
				}
				else if (vertexLinks[currentVertex].y == prevVertex)
				{
					vertexLinks[currentVertex].y = -1;
					vertexLinks[currentVertex].count--;
				}
				prevVertex = currentVertex;
				currentVertex = (long)vertexLinks[currentVertex].x;
			}
			else
			{
				break;//this is an error, and shouldn't happen
			}
		}

		if (currentVertex == -1 || currentVertex >= sizeOfVertexLinks)
			continue; //this is an error

		finalVertexLinks[count1vertices[i]].v[finalVertexLinks[count1vertices[i]].count] = vertices[currentVertex]->coord;

		finalVertexLinks[count1vertices[i]].count++;

		if (vertexLinks[currentVertex].count == 1)
		{
			vertexLinks[currentVertex].count = 0;
			vertexLinks[currentVertex].x = -1;
		}
		else if (vertexLinks[currentVertex].count == 3)
		{
			if (vertexLinks[currentVertex].x == prevVertex)
			{
				vertexLinks[currentVertex].x = -1;
			}
			else if (vertexLinks[currentVertex].y == prevVertex)
			{
				vertexLinks[currentVertex].y = -1;
			}
			else if (vertexLinks[currentVertex].z == prevVertex)
			{
				vertexLinks[currentVertex].z = -1;
			}
		}
	}
	//at this stage, all the edges that end in leaf nodes are processed, now just do the edges between vertices with 3 connections
	for (i = 0; i < count3counter; i++)
	{
		//get the (possible) three vertices that the vertex at pos count3vertices[i] of vertexLinks is linked to
		count3links[0] = (long)vertexLinks[count3vertices[i]].x;
		count3links[1] = (long)vertexLinks[count3vertices[i]].y;
		count3links[2] = (long)vertexLinks[count3vertices[i]].z;

		for (j = 0; j < 3; j++)//process each of the vertices that the current one is linked to
		{
			if (count3links[j] == -1) //all links to leaf nodes are marked with -1
				continue;

			currentVertex = count3links[j];
			prevVertex = count3vertices[i];
			while (currentVertex >= 0 && currentVertex < sizeOfVertexLinks && vertexLinks[currentVertex].count != 3)
			{
				if (vertexLinks[currentVertex].count == 2)
				{
					if (vertexLinks[currentVertex].x == prevVertex)
					{
						vertexLinks[currentVertex].x = -1;
						prevVertex = currentVertex;
						currentVertex = (long)vertexLinks[currentVertex].y;
					}
					else if (vertexLinks[currentVertex].y == prevVertex)
					{
						vertexLinks[currentVertex].y = -1;
						prevVertex = currentVertex;
						currentVertex = (long)vertexLinks[currentVertex].x;
					}
					else
					{
						break;//this is an error, prevents infinite recursion in case I make a mistake
					}

				}
			}
			if (currentVertex < 0 || currentVertex >= sizeOfVertexLinks)
			{
				continue; //this is an error, and shouldn't happen
			}

			finalVertexLinks[count3vertices[i]].v[finalVertexLinks[count3vertices[i]].count] =
				vertices[currentVertex]->coord;
			finalVertexLinks[count3vertices[i]].count++;

			if (vertexLinks[currentVertex].count == 3)
			{
				if (vertexLinks[currentVertex].x == prevVertex)
				{
					vertexLinks[currentVertex].x = -1;
				}
				else if (vertexLinks[currentVertex].y == prevVertex)
				{
					vertexLinks[currentVertex].y = -1;
				}
				else if (vertexLinks[currentVertex].z == prevVertex)
				{
					vertexLinks[currentVertex].z = -1;
				}
			}
			else
			{
				continue;
			}
		}
	}
}

bool LNLib::VoronoiDiagramGenerator::getNextVertexPair(float& x1, float& y1, float& x2, float& y2)
{
	if (finalVertexLinks == 0)
		return false;

	while (currentVertexLink < sizeOfFinalVertexLinks && finalVertexLinks[currentVertexLink].count < 1)
		currentVertexLink++;

	if (currentVertexLink >= sizeOfFinalVertexLinks || finalVertexLinks[currentVertexLink].count < 1)
		return false;

	int num = finalVertexLinks[currentVertexLink].count - 1;
	finalVertexLinks[currentVertexLink].count--;

	x1 = finalVertexLinks[currentVertexLink].v[num].x;
	y1 = finalVertexLinks[currentVertexLink].v[num].y;
	x2 = finalVertexLinks[currentVertexLink].coord.x;
	y2 = finalVertexLinks[currentVertexLink].coord.y;

	return true;
}

int LNLib::scomp(const void* p1, const void* p2)
{
	struct PointVDG* s1 = (PointVDG*)p1, * s2 = (PointVDG*)p2;
	if (s1->y < s2->y) return(-1);
	if (s1->y > s2->y) return(1);
	if (s1->x < s2->x) return(-1);
	if (s1->x > s2->x) return(1);
	return(0);
}

/* return a single in-storage site */
struct LNLib::Site* LNLib::VoronoiDiagramGenerator::nextone()
{
	struct LNLib::Site* s;
	if (siteidx < nsites)
	{
		s = &sites[siteidx];
		siteidx += 1;
		return(s);
	}
	else
		return((struct LNLib::Site*)NULL);
}
