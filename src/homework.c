# include "tsunami.h"


double interpolate(double *phi, double *U, int *map, int n){

    double u = 0.0; int i;
    for (i=0; i <n; i++)
        u += phi[i]*U[map[i]];
    return u;
    
}


void stereoCoordonnee(double* X, double* Y, double* x, double* y, double* z, int nNode){

    for(int i=0; i<nNode; i++){
        X[i] = (2*R*x[i])/(R+z[i]);
        Y[i] = (2*R*y[i])/(R+z[i]);
    }

}

void femShallowTriangleMap(int index, int map[3])
{
    int j;
    for (j=0; j < 3; ++j) 
        map[j] = index*3 + j; 
}

void femShallowEdgeMap(int index, int map[2][2], femEdge edges, int nElem, int* elem)
{
    int i,j,k;   
    for (j=0; j < 2; ++j) {
        int node = edges[index].node[j];
        for (k=0; k < 2; k++) {
            int elem = edges[index].elem[k];
            map[k][j] = (nElem)*3;
            if (elem >= 0) {
                for (i=0; i < 3; i++) {
                    if (elem[elem*3 + i] == node) {
                        map[k][j] = elem*3 + i;  }}}}}
}

void phiCreate1(double xsi, double eta, double *phi)
{
    phi[0] = 1 - xsi - eta;  
    phi[1] = xsi;
    phi[2] = eta;
}

void phiCreate2(double xsi, double *phi) 
{   
    phi[0] = (1.0 - xsi)/2.0;
    phi[1] = (1.0 + xsi)/2.0;   
    
}

void femShallowAddIntegralsElements(int nElem, int* elem, double* X, double* Y, double* E, double* U, double* V, double* BE, double* BU, double* BV){

    int iElem, i, j, k, mapElem[3];
    double xsi,eta,weight,jac,xLoc[3],yLoc[3],phi[3],dphidx[3],dphidy[3];
    double y,u,v,e;

    for (iElem=0; iElem<nElem; iElem++){
        femShallowTriangleMap(iElem, mapElem);
        int *mapCoord = &(elem[iElem*3]);
        for(j=0; j<3; j++){
            xLoc[j] = X[mapCoord[j]];
            yLoc[j] = Y[mapCoord[j]];
        }
        jac = (xLoc[1] - xLoc[0]) * (yLoc[2] - yLoc[0]) - (yLoc[1] - yLoc[0]) * (xLoc[2] - xLoc[0]);
        dphidx[0] = (yLoc[1] - yLoc[2])/jac;
        dphidx[1] = (yLoc[2] - yLoc[0])/jac;
        dphidx[2] = (yLoc[0] - yLoc[1])/jac;
        dphidy[0] = (xLoc[2] - xLoc[1])/jac;
        dphidy[1] = (xLoc[0] - xLoc[2])/jac;
        dphidy[2] = (xLoc[1] - xLoc[0])/jac;
        for (k=0; k<3; k++){
            xsi = gaussTriangleXsi[k];
            eta = gaussTriangleEta[k];
            weight = gaussTriangleWeight[k];
            phiCreate1(xsi, eta, phi);
            y = interpolate(phi, Y, mapCoord, 3);
            e = interpolate(phi, E, mapElem, 3);
            u = interpolate(phi, U, mapElem, 3); //Est-ce qu'on a déjà U et V ?
            v = interpolate(phi, V, mapElem, 3);
            double coriolis; //Je vois pas trop comment calculer ces termes
            double tau; //pareil
            double h; //Bathymétrie ??
            double toAdd = (4*R*R+xLoc[k]*xLoc[k]+yLoc[k]*yLoc[k])/(4*R*R); //A vérifier avec xLoc[k]?
            for(i=0; i<3; i++){
                BE[mapElem[i]] += ( h*u*dphidx[i] + h*v*dphidy[i] )*jac*weight; //Il faut encore ajouter les termes du tsunami
                BU[mapElem[i]] += ((coriolis*v + tau - Gamma*u)*phi[i] + g*e*dphidx[i])*jac*weight; // Pour l'instant copier coller du devoir 9
                BV[mapElem[i]] += ((- coriolis*u - Gamma*v)*phi[i] + g*e*dphidy[i])*jac*weight;}}} //Idem
            }


        }

    }
            
}

void femShallowAddIntegralsEdges(int nEdges, femEdge* edges, int nElem, int* elem, double* X, double* Y, double* E, double* U, double* V, double* BE, double* BU, double* BV){

    int iEdge, mapEdge[2][2], j, k;
    double xEdge[2], yEdge[2], phiEdge[2];
    double xsi,weight,jac;

    for(iEdge=0; iEdge<nEdges; iEdge++){
        femShallowEdgeMap(iEdge, mapEdge, edges, nElem, elem);
        for (j=0; j<2; j++){
            int node = edges[iEdge].node[j];
            xEdge[j] = X[node];
            yEdge[j] = Y[node];
        }

        //On trouve où sizeGlo ??
        int boundary = (mapEdge[1][0] == sizeGlo-1);

        double dxdxsi = (xEdge[1] - xEdge[0]);
        double dydxsi = (yEdge[1] - yEdge[0]);
        double norm = sqrt(dxdxsi*dxdxsi + dydxsi*dydxsi);
        double nx =  dydxsi/norm;
        double ny = -dxdxsi/norm;
        jac = norm / 2.0;
        for(k=0; k<2; k++){
            xsi = gaussEdgeXsi[k];
            weight = gaussEdgeWeight[k];
            phiCreate2(xsi, phiEdge);
                        
            //A vérifier
            eL = interpolate(phiEdge,E,mapEdge[0],2);
            eR = boundary ? eL : interpolate(phiEdge,E,mapEdge[1],2);
            uL = interpolate(phiEdge,U,mapEdge[0],2);
            uR = interpolate(phiEdge,U,mapEdge[1],2);
            vL = interpolate(phiEdge,V,mapEdge[0],2);
            vR = interpolate(phiEdge,V,mapEdge[1],2);
            unL = uL*nx+ vL*ny;
            unR = boundary ? -unL : uR*nx + vR*ny;
            qe =  0.5*h*   ( (unL+unR) + sqrt(g/h)*( eL-eR ) );
            qu =  0.5*g*nx*( ( eL+eR ) + sqrt(h/g)*(unL-unR) );
            qv =  0.5*g*ny*( ( eL+eR ) + sqrt(h/g)*(unL-unR) );        
            for (i=0; i < 2; i++) {
                //Termes à rajouter ?
                BE[mapEdge[0][i]] -= qe*phiEdge[i]*jac*weight; 
                BU[mapEdge[0][i]] -= qu*phiEdge[i]*jac*weight; 
                BV[mapEdge[0][i]] -= qv*phiEdge[i]*jac*weight; 
                BE[mapEdge[1][i]] += qe*phiEdge[i]*jac*weight;
                BU[mapEdge[1][i]] += qu*phiEdge[i]*jac*weight;
                BV[mapEdge[1][i]] += qv*phiEdge[i]*jac*weight; }}}

        }

    }

}




void tsunamiCompute(double dt, int nmax, int sub, const char *meshFileName, const char *baseResultName)
{ 

//
// A modifer :-)
// Ici, on se contente d'initialiser le champs d'�l�vation avec la bathym�trie
// pour avoir un joli plot de d�part !
//
// On divise la bathym�trie par un facteur (10*BathMax) pour que l'�chelle de couleur fonctionne
// bien dans la fonction d'animation fournie....
//

    int i,j,nNode,nElem, nEdges, trash, iElem, mapElem[3];
    double dtrash, xLoc[3], yLoc[3], phi[3], dphidx[3], dphidy[3];
    
    double BathMax = 9368;
    
    FILE* file = fopen(meshFileName,"r");
    fscanf(file, "Number of nodes %d \n",&nNode);   
    double *bath = malloc(sizeof(double)*nNode);
    double *x = malloc(sizeof(double)*nNode);
    double *y = malloc(sizeof(double)*nNode);
    for (i = 0; i < nNode; i++) 
        fscanf(file,"%d : %le %le %le\n",&trash,&x[i],&y[i],&bath[i]); 
    fscanf(file, "Number of triangles %d \n",&nElem); 
    int *elem = malloc(sizeof(int)*3*nElem);
    for (i = 0; i < nElem; i++) 
        fscanf(file,"%d : %d %d %d \n",&trash,&elem[i*3],&elem[i*3+1],&elem[i*3+2]);   
    fscanf(file, "Number of edges %d \n", &nEdges);
    femEdge *edges = malloc(sizeof(femEdge)*nEdges);
    for (i = 0; i<nEdges; i++){
        fscanf(file, "%d : %d %d : %d %d \n", &trash, &(edges[i].node[0]), &(edges[i].node[1]), &(edges[i].elem[0]), &(edges[i].elem[1]));
        fprintf(stdout, "%d %d \n", &(edges[i].node[0]), &(edges[i].node[1]));
    }
    double *X = malloc(sizeof(double)*nNode);
    double *Y = malloc(sizeof(double)*nNode);
    stereoCoordonnee(X, Y, x, y, bath);


    fclose(file); 


    
    double *E  = malloc(sizeof(double)*nElem*3);
    for (i = 0; i < nElem; i++)
        for (j = 0; j < 3; j++)
            E[i*3+j] = bath[elem[i*3+j]]/(10*BathMax);


    //appeler les 2 fonctions femShallowAddIntegralsElements et femShallowAddIntegralsEdges
  
    tsunamiWriteFile(baseResultName,0,E,E,E,nElem,3); 
    
    free(bath);
    free(E);
    free(elem);
 
}








