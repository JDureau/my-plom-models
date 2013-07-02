/**************************************************************************
 *    This file is part of plom.
 *
 *    plom is free software: you can redistribute it and/or modify it
 *    under the terms of the GNU General Public License as published
 *    by the Free Software Foundation, either version 3 of the
 *    License, or (at your option) any later version.
 *
 *    plom is distributed in the hope that it will be useful, but
 *    WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public
 *    License along with plom.  If not, see
 *    <http://www.gnu.org/licenses/>.
 *************************************************************************/

#include "simulation.h"

/* automatically generated code: order of the parameters */

#define ORDER_S 0
#define ORDER_I 1
#define ORDER_r0 2
#define ORDER_v 3
#define ORDER_sto 4
#define ORDER_i 5
#define ORDER_rep 6
#define ORDER_phi 7


#define ORDER_U 2
#define ORDER_R 3




#define ORDER_N 0
#define ORDER_mu_b 1
#define ORDER_mu_d 2
#define ORDER_prop 3

int step_lyap (double t, const double X[], double f[], void *params)
{
    struct s_calc *p_calc = (struct s_calc *) params;
    struct s_par *p_par = p_calc->p_par;  /* syntaxic shortcut */
    struct s_data *p_data = p_calc->p_data;

    int i;
    int c, ac, cac;    
    double **par = p_par->natural;

    struct s_router **routers = p_data->routers;  /* syntaxic shortcut */

    /* non linear system (automatically generated code)*/
    double _r[N_CAC][4];

    

    

    for(cac=0;cac<N_CAC;cac++) {
	

        

        
        _r[cac][0] = par[ORDER_r0][routers[ORDER_r0]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(X[ORDER_I*N_CAC+cac]+par[ORDER_i][routers[ORDER_i]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][1] = gsl_spline_eval(p_calc->spline[ORDER_mu_b][cac],t,p_calc->acc[ORDER_mu_b][cac])*gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][2] = gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac]);
        _r[cac][3] = (par[ORDER_v][routers[ORDER_v]->map[cac]]);
    }

    for (c=0;c<N_C;c++) {
        for(ac=0; ac<N_AC; ac++) {
            cac = c*N_AC+ac;

	    
	    f[0*N_CAC+cac] =  - (_r[cac][0]*X[ORDER_S*N_CAC+cac]) - (_r[cac][2]*X[ORDER_S*N_CAC+cac]) + (_r[cac][1]);
	    f[1*N_CAC+cac] =  - (_r[cac][3]*X[ORDER_I*N_CAC+cac]) - (_r[cac][2]*X[ORDER_I*N_CAC+cac]) + (_r[cac][0]*X[ORDER_S*N_CAC+cac]);
        }
    }

    /* linear system: product of jacobian matrix (DIM*DIM) per

       | y[1*DIM+0]       y[1*DIM+1] ...     y[1*DIM+(DIM-1)]   |
       | y[2*DIM+0]       y[2*DIM+1] ...     y[2*DIM+(DIM-1)]   |
       | ...                                                    |
       | y[DIM*DIM+0]     y[DIM*DIM+1] ...  y[DIM*DIM+(DIM-1)]  |

       (automaticaly generated code)
    */

    double _rj[N_CAC][4];
    for(cac=0; cac<N_CAC; cac++){
        
        _rj[cac][0] = -gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac])-par[ORDER_r0][routers[ORDER_r0]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(X[ORDER_I*N_CAC+cac]+par[ORDER_i][routers[ORDER_i]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][1] = -gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac])-(par[ORDER_v][routers[ORDER_v]->map[cac]])+X[ORDER_S*N_CAC+cac]*par[ORDER_r0][routers[ORDER_r0]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][2] = -X[ORDER_S*N_CAC+cac]*par[ORDER_r0][routers[ORDER_r0]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][3] = par[ORDER_r0][routers[ORDER_r0]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(X[ORDER_I*N_CAC+cac]+par[ORDER_i][routers[ORDER_i]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
    }


    for(c=0; c<N_C; c++) {
        for(ac=0; ac<N_AC; ac++) {
            cac = c*N_AC+ac;
            for(i=0; i<(N_PAR_SV*N_CAC); i++) {
                
                //printf("%d %d %d %d\n", N_PAR_SV*N_CAC, (0*N_CAC+ cac)*N_PAR_SV*N_CAC, i,  N_PAR_SV*N_CAC+ (0*N_CAC+ cac)*N_PAR_SV*N_CAC +i);
                f[N_PAR_SV*N_CAC+ (0*N_CAC+ cac)*N_PAR_SV*N_CAC +i] = +(_rj[cac][0])*X[N_PAR_SV*N_CAC+ (0*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][2])*X[N_PAR_SV*N_CAC+ (1*N_CAC+ cac)*N_PAR_SV*N_CAC +i];
                
                //printf("%d %d %d %d\n", N_PAR_SV*N_CAC, (1*N_CAC+ cac)*N_PAR_SV*N_CAC, i,  N_PAR_SV*N_CAC+ (1*N_CAC+ cac)*N_PAR_SV*N_CAC +i);
                f[N_PAR_SV*N_CAC+ (1*N_CAC+ cac)*N_PAR_SV*N_CAC +i] = +(_rj[cac][3])*X[N_PAR_SV*N_CAC+ (0*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][1])*X[N_PAR_SV*N_CAC+ (1*N_CAC+ cac)*N_PAR_SV*N_CAC +i];
                
            }
        }
    }

    return GSL_SUCCESS;
}
