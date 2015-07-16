#include "headers.h"
#include "init_2.h"
#include "declarations_2.h"

/************************************************************************/
//First derivative in X-direction 

void tdma1x(fval * fvar,int ind_var)
{	
	int ind; 

	vector<double> Qx(nx+1), Qstarx(nx+1);
	vector<double> phi(nx+1), dphidx(nx+1);
	vector<double> MWx(nx+1), MEx(nx+1), MPx(nx+1);	

	for(int j=0;j<=ny;j++)
	{
		//Coefficients for the x-direction 
		for(int i=0;i<=nx;i++)
		{
			if(i==0)
			{
				MPx[i] = 1.0; 
				MEx[i] = 2.0;
				MWx[i] = 0.0; 
			}
			else if(i==1)
			{
				MPx[i] = 1.0; 
				MWx[i] = 1./4.; 
				MEx[i] = 1./4.; 
			}
			else if(i==nx-1)
			{
				MPx[i] = 1.0; 
				MWx[i] = 1./4.; 
				MEx[i] = 1./4.; 
			}
			else if(i==nx)
			{
				MPx[i] = 1.0; 
				MWx[i] = 2.0; 
				MEx[i] = 0.0;
			}
			else
			{
				MWx[i]=1./3.; 
				MEx[i]=1./3.; 
				MPx[i]=1.; 
			}
		}
		

		for(int i=0;i<=nx;i++)
		{	
			ind = i + j*str_x; 
	
			phi[i]=fvar[ind].u[ind_var];
			dphidx[i]=0.0; 
			Qx[i] = 0.0; 
			Qstarx[i] = 0.0; 
		}

		//Populating the matrix

		for(int i=0;i<=nx;i++)
		{

			if(i==0)
			{
				Qx[i] = (1./dx)*( (-15./6.)*phi[i] + 2.*phi[i+1] + (1./2.)*phi[i+2]);
			}
			else if(i==1)
			{				
				Qx[i] = (3./2.)*( phi[i+1] - phi[i-1])/(2.*dx); 
			}
			else if(i==nx-1)
			{				
				Qx[i] = (3./2.)*( phi[i+1] - phi[i-1])/(2.*dx); 
			}
			else if(i==nx)
			{ 
				Qx[i] = (1./dx)*( (15./6.)*phi[i] - 2.*phi[i-1] - (1./2.)*phi[i-2])  ; 
			}
			else
			{
				Qx[i] = (phi[i+2] - phi[i-2])*(td_bx/(4.*dx)) + (phi[i+1]-phi[i-1])*(td_ax/(2.*dx)); 
			}

		}

// Thomas Algorithm for X-derivative 

//Forward Elimination

		for(int i=0;i<=nx;i++)
		{
			if(i==0)
			{
				Qstarx[i]=Qx[i];
			}
			else 
			{
				MPx[i] =  MPx[i]- ( MWx[i]*MEx[i-1] )/(MPx[i-1]) ; 
				Qstarx[i] = Qx[i] - ( MWx[i]*Qstarx[i-1] )/(MPx[i-1]) ; 
			}
		}

//Backward Substitution
	
		for(int i=nx;i>=0;i--)
		{	
			if(i==nx)
			{
				dphidx[i] = Qstarx[i]/MPx[i]; 
			}	
			else 
			{
				dphidx[i] = ( Qstarx[i] - MEx[i]*dphidx[i+1] )/MPx[i]; 
			}
	
			ind = i + j*str_x; 
		
			fvar[ind].ux[ind_var] = dphidx[i];		
		}

	} 

}

/************************************************************************/
//First derivative in Y-direction 

void tdma1y(fval * fvar, int ind_var)
{
	int ind; 

	vector<double> Qy(ny+1), Qstary(ny+1);
	vector<double> phi(ny+1), dphidy(ny+1);
	vector<double> MWy(ny+1), MEy(ny+1), MPy(ny+1);	

	for(int i=0;i<=nx;i++)
	{
		//Coefficients for the y-direction 
		for(int j=0;j<=ny;j++)
		{
			if(j==0)
			{
				MPy[j] = 1.0; 
				MEy[j] = 2.0;
				MWy[j] = 0.0; 
			}
			else if(j==1)
			{
				MPy[j] = 1.0; 
				MWy[j] = 1./4.; 
				MEy[j] = 1./4.; 
			}
			else if(j==ny-1)
			{
				MPy[j] = 1.0; 
				MWy[j] = 1./4.; 
				MEy[j] = 1./4.; 
			}
			else if(j==ny)
			{
				MPy[j] = 1.0; 
				MWy[j] = 2.0; 
				MEy[j] = 0.0;
			}
			else
			{
				MWy[j]=1./3.; 
				MEy[j]=1./3.; 
				MPy[j]=1.; 
			}
		}

		for(int j=0;j<=ny;j++)
		{	
			ind = i + j*str_x; 
	
			phi[j]=fvar[ind].u[ind_var];
			dphidy[j]=0.0; 
			Qy[j] = 0.0; 
			Qstary[j] = 0.0; 
		}

		//*Populating the matrix

		for(int j=0;j<=ny;j++)
		{

			if(j==0)
			{
				Qy[j] = (1./dy)*( (-15./6.)*phi[j] + 2.*phi[j+1] + (1./2.)*phi[j+2] );
			}
			else if(j==1)
			{				
				Qy[j] = (3./2.)*( phi[j+1] - phi[j-1])/(2.*dy); 
			}
			else if(j==ny-1)
			{				
				Qy[j] = (3./2.)*( phi[j+1] - phi[j-1])/(2.*dy); 
			}
			else if(j==ny)
			{ 
				Qy[j] = (1./dy)*( (15./6.)*phi[j] - 2.*phi[j-1] - (1./2.)*phi[j-2] )  ; 
			}
			else
			{
				Qy[j] = (phi[j+2] - phi[j-2])*(td_by/(4.*dy)) + (phi[j+1]-phi[j-1])*(td_ay/(2.*dy)); 
			}

		}

// Thomas Algorithm for X-derivative 

//Forward Elimination

		for(int j=0;j<=ny;j++)
		{
			if(j==0)
			{
				Qstary[j]=Qy[j];
			}
			else 
			{
				MPy[j] =  MPy[j]- ( MWy[j]*MEy[j-1])/(MPy[j-1]) ; 
				Qstary[j] = Qy[j] - ( MWy[j]*Qstary[j-1] )/(MPy[j-1]) ; 
			}
		}

//Backward Substitution
	
		for(int j=ny;j>=0;j--)
		{	
			if(j==ny)
			{
				dphidy[j] = Qstary[j]/MPy[j]; 
			}	
			else 
			{
				dphidy[j] = ( Qstary[j] - MEy[j]*dphidy[j+1] )/MPy[j]; 
			}
		
			ind = i + j*str_x; 
		
			fvar[ind].uy[ind_var] = dphidy[j];		
		}

	} 
}

/************************************************************************/
//Second derivative in X-direction 

void tdma2x(fval * fvar, int ind_var)
{
	int ind; 

	vector<double> Qx(nx+1), Qstarx(nx+1);
	vector<double> phi(nx+1), dphidx(nx+1);
	vector<double> MW2x(nx+1), ME2x(nx+1), MP2x(nx+1); 
	

	//Initializing the array Phi

	for(int j=0;j<=ny;j++) 
	{ 
		//Coefficients in the x-direction 

		for(int i=0;i<=nx;i++)
		{
			if(i==0)
			{
				MP2x[i] = 1.0; 
				ME2x[i] = 11.0;
				MW2x[i] = 0.0; 
			}
			else if(i==1)
			{
				MP2x[i] = 1.0; 
				ME2x[i] = 1./10.;
				MW2x[i] = 1./10.;   
			} 
			else if(i==nx-1)
			{
				MP2x[i] = 1.0; 
				ME2x[i] = 1./10.;
				MW2x[i] = 1./10.;   
			} 
			else if(i==nx)
			{
				MP2x[i] = 1.0; 
				MW2x[i] = 11.0; 
				ME2x[i] = 0.0;
			}
			else
			{
				MW2x[i]=2./11.; 
				ME2x[i]=2./11.; 
				MP2x[i]=1.; 
			}

		}


		for(int i=0;i<=nx;i++)
		{
			ind = i + j*str_x; 
     
			phi[i]=fvar[ind].u[ind_var];
			dphidx[i]=0.0; 
			Qx[i] = 0.0; 
			Qstarx[i] = 0.0; 
		}

		//Populating the matrix

		for(int i=0;i<=nx;i++)
		{
			if(i==0)
			{
				Qx[i] = ( 1./(dx*dx) )*( 13.*phi[i] - 27.*phi[i+1] + 15.*phi[i+2] - phi[i+3]);
			}
			else if(i==1)
			{
				Qx[i] = (12./(10.*dx*dx))*( phi[i+1] - 2.*phi[i] + phi[i-1]) ;
			}
			else if(i==nx-1)
			{
				Qx[i] = (12./(10.*dx*dx))*( phi[i+1] - 2.*phi[i] + phi[i-1]) ;
			}
			else if(i==nx)
			{        
				Qx[i] = ( 1./(dx*dx) )*( 13.*phi[i] - 27.*phi[i-1] + 15.*phi[i-2] - phi[i-3]);
			}
			else
			{
				Qx[i] = (1.0/(dx*dx))*( td_ax2*(phi[i+1] - 2.*phi[i] + phi[i-1]) + (td_bx2/4.)*( phi[i+2] - 2.*phi[i] + phi[i-2] ) ) ; 
			}
		}

// Thomas Algorithm for X-derivative 

//Forward Elimination

		for(int i=0;i<=nx;i++)
		{
			if(i==0)
			{
				Qstarx[i]=Qx[i];
			}
			else 
			{
				MP2x[i] =  MP2x[i]- ( MW2x[i]*ME2x[i-1])/(MP2x[i-1]) ; 
				Qstarx[i] = Qx[i] - ( MW2x[i]*Qstarx[i-1] )/(MP2x[i-1]) ; 
			}
		}

//Backward Substitution
	
		for(int i=nx;i>=0;i--)
		{	
			if(i==nx)
			{
				dphidx[i] = Qstarx[i]/MP2x[i]; 
			}	
			else 
			{
				dphidx[i] = ( Qstarx[i] - ME2x[i]*dphidx[i+1] )/MP2x[i]; 
			}
	   
			ind = i + j*str_x; 
		
			fvar[ind].uxx[ind_var] = dphidx[i];		
		}

	} 

}

/************************************************************************/
//Second derivative in Y-direction 

void tdma2y(fval * fvar, int ind_var)
{
	double td_ay = 12./11.; 
	double td_by = 3./11.;

	int ind; 

	vector<double> Qy(ny+1), Qstary(ny+1);
	vector<double> phi(ny+1), dphidy(ny+1);
	vector<double> MW2y(ny+1), ME2y(ny+1), MP2y(ny+1); 

	//Initializing the array Phi

	for(int i=0;i<=nx;i++) 
	{ 

		//Coefficients in the y-direction 

		for(int j=0;j<=ny;j++)
		{
			if(j==0)
			{
				MP2y[j] = 1.0; 
				ME2y[j] = 11.0;
				MW2y[j] = 0.0; 
			}
			else if(j==1)
			{
				MP2y[j] = 1.0; 
				ME2y[j] = 1./10.;
				MW2y[j] = 1./10.;   
			} 
			else if(j==ny-1)
			{
				MP2y[j] = 1.0; 
				ME2y[j] = 1./10.;
				MW2y[j] = 1./10.;   
			} 
			else if(j==ny)
			{
				MP2y[j] = 1.0; 
				MW2y[j] = 11.0; 
				ME2y[j] = 0.0;
			}
			else
			{
				MW2y[j]=2./11.; 
				ME2y[j]=2./11.; 
				MP2y[j]=1.; 
			}

		}


		for(int j=0;j<=ny;j++)
		{
			ind = i + j*str_x; 
     
			phi[j]=fvar[ind].u[ind_var];
			dphidy[j]=0.0; 
			Qy[j] = 0.0; 
			Qstary[j] = 0.0; 
		}

		//Populating the matrix

		for(int j=0;j<=ny;j++)
		{

			if(j==0)
			{
				Qy[j] = (1./(dy*dy) )*( 13.*phi[j] - 27.*phi[j+1] + 15.*phi[j+2] - phi[j+3]);
			}
			else if(j==1)
			{
				Qy[j] = (12./(10.*dy*dy))*( phi[j+1] - 2.*phi[j] + phi[j-1]) ;
			}
			else if(j==ny-1)
			{
				Qy[j] = (12./(10.*dy*dy))*( phi[j+1] - 2.*phi[j] + phi[j-1]) ;
			}
			else if(j==ny)
			{        
				Qy[j] = ( 1./(dy*dy) )*( 13.*phi[j] - 27.*phi[j-1] + 15.*phi[j-2] - phi[j-3]);
			}
			else
			{
				Qy[j] = ( 1.0/(dy*dy) )*( td_ay*( phi[j+1] - 2.*phi[j] + phi[j-1] ) + (td_by/4.)*( phi[j+2] - 2.*phi[j] + phi[j-2] ) ) ; 
			}

		}

// Thomas Algorithm for Y-derivative 

//Forward Elimination

		for(int j=0;j<=ny;j++)
		{
			if(j==0)
			{
				Qstary[j]=Qy[j];
			}
			else 
			{
				MP2y[j] =  MP2y[j]- ( MW2y[j]*ME2y[j-1])/(MP2y[j-1]) ; 
				Qstary[j] = Qy[j] - ( MW2y[j]*Qstary[j-1] )/(MP2y[j-1]) ; 
			}
		}

//Backward Substitution
	
		for(int j=ny;j>=0;j--)
		{	
			if(j==ny)
			{
				dphidy[j] = Qstary[j]/MP2y[j]; 
			}	
			else 
			{
				dphidy[j] = ( Qstary[j] - ME2y[j]*dphidy[j+1] )/MP2y[j]; 
			}
	   
			ind = i + j*str_x; 
		
			fvar[ind].uyy[ind_var] = dphidy[j];		
		}

	} 

}
/***************************************************************************************************************************/
//Interpolation functions that interpolates the advection terms from the corners of cell to the edge centers of the adjacent cells in the X- and Y-directions. This is used to compute the source term of the Poisson's equation
/***************************************************************************************************************************/
//This is to interpolate from a cell corner onto the adjacent edge centers     ---X----o----X----

void interpol_cen_hface(fval* fvar, fval*  hvar, int ind_var) //Checked for correctness
{  
	int ind; 
	
	vector<double> Qx(nx), Qstarx(nx);
	vector<double> phi(nx+1), phi_polated(nx);
	vector<double> MWix(nx), MEix(nx), MPix(nx); 	

	/*******Initializing the array Phi*/
	
	for(int j=0;j<=ny;j++) 
	{ 
		//Coefficients in the x-direction 

		for(int i=0;i<nx;i++)
		{
			if(i==0)
			{								
				MPix[i]=1.0;   //4th order scheme from Knikker 
				MEix[i]=1.0;
				MWix[i]=0.0; 
			}
			else if(i==nx-1)			
			{								
				MPix[i]=1.0; 
				MWix[i]=1.0; 
				MEix[i]=0.0; 
			}
			else
			{				
				MPix[i]=6.0; 
				MEix[i]=1.0;
				MWix[i]=1.0; 				
			}
		}

		for(int i=0;i<=nx;i++)
		{
		  ind = i + j*str_x; 
                  phi[i]=fvar[ind].u[ind_var];
		}
		
		for(int i=0;i<nx;i++)
		{
		    phi_polated[i] = 0.0;
		    Qx[i] = 0.0; 
		    Qstarx[i] = 0.0; 		   
		}

		/*******Populating the matrix*/

		for(int i=0;i<nx;i++)
		{
			if(i==0)
			{								 				
				Qx[i] = (1./4.)*phi[i] + (3./2.)*phi[i+1] + (1./4.)*phi[i+2]; 
			}
			else if(i==nx-1)
			{				
				
				Qx[i] = (1./4.)*phi[i+1] + (3./2.)*phi[i] + (1./4.)*phi[i-1]; 
			}
			else
			{				
				Qx[i] = 4.0*( phi[i+1] + phi[i] ); 
			}
		}

//Forward Elimination

		for(int i=0;i<nx;i++)
		{
			if(i==0)
			{
				Qstarx[i]=Qx[i];
			}
			else 
			{
				MPix[i] =  MPix[i]- ( MWix[i]*MEix[i-1] )/(MPix[i-1]) ; 
				Qstarx[i] = Qx[i] - ( MWix[i]*Qstarx[i-1] )/(MPix[i-1]) ; 
			}
		}

//Backward Substitution
	
		for(int i=nx-1;i>=0;i--)
		{	
			if(i==nx-1)
			{
				phi_polated[i] = Qstarx[i]/MPix[i]; 
			}	
			else 
			{
				phi_polated[i] = ( Qstarx[i] - MEix[i]*phi_polated[i+1] )/MPix[i]; 
			}
	   
			ind = i+j*cstr_x; 
		
			hvar[ind].u[ind_var] = phi_polated[i];	
			
			//cout<<j<<"	"<<i<<"		"<<hvar[ind].u[ind_var]<<"\n";	
		}
	}
}

/***************************************************************************************************************************/
//This is to interpolate from a cell corner onto the adjacent vertical edge centers
//	 0
//	 |
//	 X
//	 |
//-o--X--o--X--o-
//	 |
//	 X
//	 |
//	 0
void interpol_cen_vface(fval* fvar, fval* vvar, int ind_var) //Checked for correctness
{  
	int ind; 
	
	vector<double> Qy(ny), Qstary(ny);
	vector<double> phi(ny+1), phi_polated(ny);
	vector<double> MWjy(ny), MEjy(ny), MPjy(ny); 	

	/*******Initializing the array Phi*/
	
	for(int i=0;i<=nx;i++) 
	{ 
		//Coefficients in the x-direction 

		for(int j=0;j<ny;j++)
		{
			if(j==0)
			{								
				MPjy[j]=1.0;   //4th order scheme from Knikker 
				MEjy[j]=1.0;
				MWjy[j]=0.0; 
			}
			else if(j==ny-1)			
			{								
				MPjy[j]=1.0; 
				MWjy[j]=1.0; 
				MEjy[j]=0.0; 
			}
			else
			{				
				MPjy[j]=6.0; 
				MEjy[j]=1.0;
				MWjy[j]=1.0; 				
			}
		}

		for(int j=0;j<=ny;j++)
		{
		  ind = i + j*str_x; 
                  phi[j]=fvar[ind].u[ind_var];
		}
		
		for(int j=0;j<ny;j++)
		{
		    phi_polated[j] = 0.0;
		    Qy[j] = 0.0; 
		    Qstary[j] = 0.0; 		   
		}

		/*******Populating the matrix*/

		for(int j=0;j<ny;j++)
		{
			if(j==0)
			{								 				
				Qy[j] = (1./4.)*phi[j] + (3./2.)*phi[j+1] + (1./4.)*phi[j+2]; 
			}
			else if(j==ny-1)
			{				
				
				Qy[j] = (1./4.)*phi[j+1] + (3./2.)*phi[j] + (1./4.)*phi[j-1]; 
			}
			else
			{				
				Qy[j] = 4.0*( phi[j+1] + phi[j] ); 
			}
		}

//Forward Elimination

		for(int j=0;j<ny;j++)
		{
			if(j==0)
			{
				Qstary[j]=Qy[j];
			}
			else 
			{
				MPjy[j] =  MPjy[j]- ( MWjy[j]*MEjy[j-1] )/(MPjy[j-1]) ; 
				Qstary[j] = Qy[j] - ( MWjy[j]*Qstary[j-1] )/(MPjy[j-1]) ; 
			}
		}

//Backward Substitution
	
		for(int j=ny-1;j>=0;j--)
		{	
			if(j==ny-1)
			{
				phi_polated[j] = Qstary[j]/MPjy[j]; 
			}	
			else 
			{
				phi_polated[j] = ( Qstary[j] - MEjy[j]*phi_polated[j+1] )/MPjy[j]; 
			}
	   
			ind = i+j*str_x; 
		
			vvar[ind].u[ind_var] = phi_polated[j];		
			
			//cout<<i<<"	"<<j<<"		"<<vvar[ind].u[ind_var]<<"\n";	
		}
	}
}

/******************************************************************************************************/
/******************Derivatives on the cell corners required for the iteration**************************/ 
/******************************************************************************************************/
//X-derivative on the corner based on the values at the adjacent edge center ---o----X----o----

void derv_hface2node_x(fval* fvar, fval* hvar, int ind_var)
{
	int ind; 
	
	vector<double> Qx(nx-1), Qstarx(nx-1);
	vector<double> phi(nx), phi_polated(nx-1);
	vector<double> MWix(nx-1), MEix(nx-1), MPix(nx-1); 	

	/*******Initializing the array Phi*/
	
	for(int j=0;j<=ny;j++) 
	{ 
		//Coefficients in the x-direction 
		for(int i=0;i<nx-1;i++)
		{
			if(i==0)
			{
			   MPix[i] = 1.0; 
			   MEix[i] = -1.0;
			   MWix[i] = 0.0; 
			}
			else if(i==nx-2)
			{
			    MPix[i] = 1.0;
			    MWix[i] = -1.0;
			    MEix[i] = 0.0;		    			
			}
			else
			{
			    MPix[i] = 1.0; 
			    MEix[i] = 9./62.;
			    MWix[i] = 9./62.;
			    
			    /*MPix[i] = 22.0; 
			    MEix[i] = 1.0;
			    MWix[i] = 1.0;*/
			}
		}

		for(int i=0;i<nx;i++)
		{
		  ind = i + j*cstr_x; 
                  phi[i] = hvar[ind].u[ind_var];
		}
		
		for(int i=0;i<nx-1;i++)
		{
		    phi_polated[i] = 0.0; 
		    Qx[i] = 0.0; 
		    Qstarx[i] = 0.0; 		   
		}

		/*******Populating the matrix*/		
		for(int i=0;i<nx-1;i++)
		{
			if(i==0)
			{
			    Qx[i] = (1./dx)*( -phi[i] + 2.*phi[i+1] - phi[i+2] );  			    
			}			
			else if(i==nx-2)
			{			
			    Qx[i] = -(1./dx)*( -phi[i+1] + 2.*phi[i] - phi[i-1] );  			
			}
			else
			{
			    Qx[i] = (1./dx)*( (63./62.)*(phi[i+1]-phi[i])  + (1./3.)*(17./62.)*(phi[i+2]-phi[i-1]) ) ;
			    
			   // Qx[i] = (24./dx)*( phi[i+1]-phi[i] ); //Taken from Knikker. 4th order in the middle and 3rd order at the end.  
			}
		}

// Thomas Algorithm for X-derivative 

//Forward Elimination

		for(int i=0;i<nx-1;i++)
		{
			if(i==0)
			{
				Qstarx[i]=Qx[i];
			}
			else 
			{
				MPix[i] =  MPix[i]- ( MWix[i]*MEix[i-1] )/(MPix[i-1]) ; 
				Qstarx[i] = Qx[i] - ( MWix[i]*Qstarx[i-1] )/(MPix[i-1]) ; 
			}
		}

//Backward Substitution
	
		for(int i=nx-2;i>=0;i--)
		{	
			if(i==nx-2)
			{
				phi_polated[i] = Qstarx[i]/MPix[i]; 
			}	
			else 
			{
				phi_polated[i] = ( Qstarx[i] - MEix[i]*phi_polated[i+1] )/MPix[i]; 
			}
	   
			ind = i+1+j*str_x; 
		
			fvar[ind].ux[ind_var] = phi_polated[i];		
			
			//cout<<j<<"	"<<i<<"		"<<fvar[ind].ux[ind_var]<<"\n";		
		}
		
		//cout<<"-------------------------------------------------------------------\n";
	}
}

/**************************************************************************/
//Y-derivative on the corner based on the values at the adjacent edge centers
//	 0
//	 |
//	 X
//	 |
//-o--X--o--X--o-
//	 |
//	 X
//	 |
//	 0	
void derv_vface2node_y(fval* fvar, fval* vvar, int ind_var)  //d/dy(vface_value) on cell center 
{
	int ind; 
	
	vector<double> Qy(ny-1), Qstary(ny-1);
	vector<double> phi(ny), phi_polated(ny-1);
	vector<double> MWy(ny-1), MEy(ny-1), MPy(ny-1); 	

	/*******Initializing the array Phi*/
	
	for(int i=1;i<=nx;i++) 
	{ 
		
		for(int j=0;j<ny-1;j++)
		{
			if(j==0)
			{
			   MPy[j] = 1.0; 
			   MEy[j] = -1.0;
			   MWy[j] = 0.0; 
			}
			else if(j==ny-2)
			{
			    MPy[j] = 1.0;
			    MWy[j] = -1.0;
			    MEy[j] = 0.0;		    			
			}
			else
			{
			    MPy[j] = 1.0; 
			    MEy[j] = 9./62.;
			    MWy[j] = 9./62.;   			  		
			}
		}


		for(int j=0;j<ny;j++) //Assigning the Vface values to the "phi" variable here
		{
		  ind =  i + j*str_x; 
                  phi[j]= vvar[ind].u[ind_var];
		}
		
		for(int j=0;j<ny-1;j++)
		{
		    phi_polated[j] = 0.0; 
		    Qy[j] = 0.0; 
		    Qstary[j] = 0.0; 		   
		}

		/*******Populating the matrix*/		
		for(int j=0;j<ny-1;j++)
		{
			if(j==0)
			{
			    Qy[j] = (1./dy)*( -phi[j] + 2.*phi[j+1] - phi[j+2] );  
			}			
			else if(j==ny-2)
			{			
			    Qy[j] = -(1./dy)*( -phi[j+1] + 2.*phi[j] - phi[j-1] );  			
			}
			else
			{
			    Qy[j] = (1./dy)*( (63./62.)*(phi[j+1]-phi[j])  + (1./3.)*(17./62.)*(phi[j+2]-phi[j-1]) ) ; 
			}
		}

// Thomas Algorithm for X-derivative 

//Forward Elimination

		for(int j=0;j<ny-1;j++)
		{
			if(j==0)
			{
				Qstary[j]=Qy[j];
			}
			else 
			{
				MPy[j] =  MPy[j]- ( MWy[j]*MEy[j-1] )/(MPy[j-1]) ; 
				Qstary[j] = Qy[j] - ( MWy[j]*Qstary[j-1] )/(MPy[j-1]) ; 
			}
		}

//Backward Substitution
	
		for(int j=ny-2;j>=0;j--)
		{	
			if(j==ny-2)
			{
				phi_polated[j] = Qstary[j]/MPy[j]; 
			}	
			else 
			{
				phi_polated[j] = ( Qstary[j] - MEy[j]*phi_polated[j+1] )/MPy[j]; 
			}
	   
			ind = i+(j+1)*cstr_x; 
		
			fvar[ind].uy[ind_var] = phi_polated[j];	
			
			//cout<<i<<"	"<<j<<"		"<<fvar[ind].uy[ind_var]<<"\n";		
		}
		
		//cout<<"------------------------------------------------------------------------\n";
	}
	

}
/*********************************************************************************/
