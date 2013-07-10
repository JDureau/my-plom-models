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

#include "plom.h"

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

/**
 * Alloc memory for the psr implementation
 */
void build_psr(struct s_calc *p)
{
    unsigned int tab[4]; 

    /*automaticaly generated code: dimension of prob and inc*/
    
    tab[ORDER_S] = 3;
    tab[ORDER_I] = 3;
    tab[ORDER_U] = 2;
    tab[ORDER_R] = 1;

    p->prob = init2d_var_set0(4, tab);
    p->inc = init3u_varp2_set0(4, N_CAC, tab);

    //  p->gravity = init1d_set0(N_C);
}


void proj2obs(struct s_X *p_X, struct s_data *p_data)
{
    int o, ind_obs, ind_proj_inc;
    int n_ts_unique_o, n_stream_o_ts;

    struct s_obs2ts **obs2ts = p_data->obs2ts;

    

    ind_obs = 0;
    ind_proj_inc = N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot;

    /* extend incidence: duplicate ts with multiple data streams */
    for(o=0; o<N_OBS_INC; o++) {
        for(n_ts_unique_o=0; n_ts_unique_o< (obs2ts[o])->n_ts_unique; n_ts_unique_o++) {
            for(n_stream_o_ts=0; n_stream_o_ts< (obs2ts[o])->n_stream[n_ts_unique_o]; n_stream_o_ts++) {
                p_X->obs[ind_obs++] = p_X->proj[ind_proj_inc];
            }
            ind_proj_inc++;
        }
    }

    /* add prevalence: aggregate across c and ac to match ts and repeat to tacle multiple data streams */
    
}

//stepping functions for Poisson System with stochastic rates (psr)
void step_psr(struct s_X *p_X, double t, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc)
{
    /* t is the time in unit of the data */

    struct s_obs2ts **obs2ts = p_data->obs2ts;  /* syntaxic shortcut */
    struct s_router **routers = p_data->routers;   /* syntaxic shortcut */

    int c, ac, cac, n_cac, ts, o;
    double sum_inc = 0.0;
    int offset;

    double sum, one_minus_exp_sum;

    double **par = p_par->natural;   

    double *X = p_X->proj;
    double dt = p_X->dt;


    /*automaticaly generated code:*/
    /*0-declaration of noise terms (if any)*/
    
    double white_noise__0;

    double _r[3];
    
    double _sf[1];

    


    for(c=0;c<N_C;c++) {
        for(ac=0;ac<N_AC;ac++) {
            cac = c*N_AC+ac;

            

            /*1-generate noise increments (if any) (automaticaly generated code)*/
            
            if(p_data->noises_off & PLOM_NO_ENV_STO){
                
                white_noise__0 = 1.0;
            } else {
                
                white_noise__0 = gsl_ran_gamma(p_calc->randgsl, (dt)/ pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]], 2), pow(par[ORDER_sto][routers[ORDER_sto]->map[cac]], 2))/dt;
            }
            

            /*2-generate process increments (automaticaly generated code)*/
            
            _sf[0] = correct_rate(par[ORDER_v][routers[ORDER_v]->map[cac]],dt);

            
            _r[0] = _sf[0];
            _r[1] = par[ORDER_r0][routers[ORDER_r0]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*white_noise__0*(X[ORDER_I*N_CAC+cac]+par[ORDER_i][routers[ORDER_i]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
            _r[2] = gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac]);

            sum = _r[1]*dt+_r[2]*dt;
if(sum>0.0){
one_minus_exp_sum = (1.0-exp(-sum));
p_calc->prob[ORDER_S][0] = one_minus_exp_sum*((_r[1]*dt)/sum);
p_calc->prob[ORDER_S][1] = one_minus_exp_sum*((_r[2]*dt)/sum);
p_calc->prob[ORDER_S][2] = 1.0 - p_calc->prob[ORDER_S][0] - p_calc->prob[ORDER_S][1];
}
else{
p_calc->prob[ORDER_S][0] = 0.0;
p_calc->prob[ORDER_S][1] = 0.0;
p_calc->prob[ORDER_S][2] = 1.0;
}

sum = _r[0]*dt+_r[2]*dt;
if(sum>0.0){
one_minus_exp_sum = (1.0-exp(-sum));
p_calc->prob[ORDER_I][0] = one_minus_exp_sum*((_r[0]*dt)/sum);
p_calc->prob[ORDER_I][1] = one_minus_exp_sum*((_r[2]*dt)/sum);
p_calc->prob[ORDER_I][2] = 1.0 - p_calc->prob[ORDER_I][0] - p_calc->prob[ORDER_I][1];
}
else{
p_calc->prob[ORDER_I][0] = 0.0;
p_calc->prob[ORDER_I][1] = 0.0;
p_calc->prob[ORDER_I][2] = 1.0;
}



            /*3-multinomial drawn (automaticaly generated code)*/
            
            plom_ran_multinomial(p_calc->randgsl, 3, (unsigned int) X[ORDER_S*N_CAC+cac], p_calc->prob[ORDER_S], p_calc->inc[ORDER_S][cac]);
            plom_ran_multinomial(p_calc->randgsl, 3, (unsigned int) X[ORDER_I*N_CAC+cac], p_calc->prob[ORDER_I], p_calc->inc[ORDER_I][cac]);

            /*4-update state variables (automaticaly generated code)*/
            //use inc to cache the Poisson draw as thew might be re-used for the incidence computation
            
            p_calc->inc[ORDER_U][cac][0] = gsl_ran_poisson(p_calc->randgsl, (gsl_spline_eval(p_calc->spline[ORDER_mu_b][cac],t,p_calc->acc[ORDER_mu_b][cac])*gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]))*dt);

            X[ORDER_S*N_CAC+cac] = p_calc->inc[ORDER_S][cac][2] + p_calc->inc[ORDER_U][cac][0];
X[ORDER_I*N_CAC+cac] = p_calc->inc[ORDER_I][cac][2] + p_calc->inc[ORDER_S][cac][0];


        }/*end for on ac*/
    } /*end for on c*/

    /*compute incidence:integral between t and t+1 (automaticaly generated code)*/

    offset = N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot;
    
    o = 0;

    for(ts=0; ts<obs2ts[o]->n_ts_unique; ts++) {
        sum_inc = 0.0;
        for(n_cac=0; n_cac<obs2ts[o]->n_cac[ts]; n_cac++) {
            c = obs2ts[o]->cac[ts][n_cac][0];
            ac = obs2ts[o]->cac[ts][n_cac][1];
            cac = c*N_AC+ac;

            sum_inc +=  + p_calc->inc[ORDER_S][cac][0];
        }
        X[offset] += sum_inc;
        offset++;
    }
    
}



//stepping functions for ODE and SDEs



int step_ode(double t, const double X[], double f[], void *params)

{

    
    struct s_calc *p_calc = (struct s_calc *) params;
    struct s_data *p_data = p_calc->p_data;
    struct s_par *p_par = p_calc->p_par;
    

    struct s_obs2ts **obs2ts = p_data->obs2ts;
    struct s_router **routers = p_data->routers;

    int i, c, ac, cac, n_cac, ts, o;
    double sum_inc = 0.0;
    int offset;
    
    double **par = p_par->natural;

    double _r[N_CAC][4];

    

    

    


    for(cac=0;cac<N_CAC;cac++){
        


        

        
        _r[cac][0] = par[ORDER_r0][routers[ORDER_r0]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(X[ORDER_I*N_CAC+cac]+par[ORDER_i][routers[ORDER_i]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][1] = gsl_spline_eval(p_calc->spline[ORDER_mu_b][cac],t,p_calc->acc[ORDER_mu_b][cac])*gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][2] = gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac]);
        _r[cac][3] = (par[ORDER_v][routers[ORDER_v]->map[cac]]);

        
    }

    for(c=0;c<N_C;c++) {
        for(ac=0; ac<N_AC; ac++) {
            cac = c*N_AC+ac;

            /*automaticaly generated code:*/
            /*ODE system*/
            
            f[0*N_CAC+cac] =  - (_r[cac][0]*X[ORDER_S*N_CAC+cac]) - (_r[cac][2]*X[ORDER_S*N_CAC+cac]) + (_r[cac][1]);
            f[1*N_CAC+cac] =  - (_r[cac][3]*X[ORDER_I*N_CAC+cac]) - (_r[cac][2]*X[ORDER_I*N_CAC+cac]) + (_r[cac][0]*X[ORDER_S*N_CAC+cac]);
        }
    }


    //TODO: drift of the diffusion
    //for(i=N_PAR_SV*N_CAC; i<(N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot); i++){
    //    f[i] = 0.0;
    //}

    /*automaticaly generated code:*/
    /*compute incidence:integral between t and t+1*/
    offset = N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot;

    
    o = 0;

    for (ts=0; ts<obs2ts[o]->n_ts_unique; ts++) {
        sum_inc = 0.0;
        for (n_cac=0; n_cac<obs2ts[o]->n_cac[ts]; n_cac++) {
            c = obs2ts[o]->cac[ts][n_cac][0];
            ac = obs2ts[o]->cac[ts][n_cac][1];
            cac = c*N_AC+ac;

            sum_inc +=  + (_r[cac][0]*X[ORDER_S*N_CAC+cac]);
        }

        f[offset] = sum_inc;
        offset++;
    }
    

    
    return GSL_SUCCESS;
    

}


void step_sde_no_env_sto(struct s_X *p_X, double t, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc)

{

    
    double *X = p_X->proj;
    double dt = p_X->dt;
    double *f = p_calc->y_pred;
    

    struct s_obs2ts **obs2ts = p_data->obs2ts;
    struct s_router **routers = p_data->routers;

    int i, c, ac, cac, n_cac, ts, o;
    double sum_inc = 0.0;
    int offset;
    
    double **par = p_par->natural;

    double _r[N_CAC][4];

    

    
    double dem_sto__0[N_CAC];
    double dem_sto__1[N_CAC];
    double dem_sto__2[N_CAC];
    double dem_sto__3[N_CAC];
    double dem_sto__4[N_CAC];

    


    for(cac=0;cac<N_CAC;cac++){
        


        

        
        _r[cac][0] = par[ORDER_r0][routers[ORDER_r0]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(X[ORDER_I*N_CAC+cac]+par[ORDER_i][routers[ORDER_i]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][1] = gsl_spline_eval(p_calc->spline[ORDER_mu_b][cac],t,p_calc->acc[ORDER_mu_b][cac])*gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][2] = gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac]);
        _r[cac][3] = (par[ORDER_v][routers[ORDER_v]->map[cac]]);

        
        dem_sto__0[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__1[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__2[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__3[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__4[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
    }

    for(c=0;c<N_C;c++) {
        for(ac=0; ac<N_AC; ac++) {
            cac = c*N_AC+ac;

            /*automaticaly generated code:*/
            /*ODE system*/
            
            f[0*N_CAC+cac] = X[0*N_CAC+cac] +  ( - (_r[cac][0]*X[ORDER_S*N_CAC+cac]) - (_r[cac][2]*X[ORDER_S*N_CAC+cac]) + (_r[cac][1]))*dt + - sqrt((_r[cac][0]*X[ORDER_S*N_CAC+cac]))*dem_sto__1[cac]- sqrt((_r[cac][2]*X[ORDER_S*N_CAC+cac]))*dem_sto__3[cac]+ sqrt((_r[cac][1]))*dem_sto__0[cac];
            f[1*N_CAC+cac] = X[1*N_CAC+cac] +  ( - (_r[cac][3]*X[ORDER_I*N_CAC+cac]) - (_r[cac][2]*X[ORDER_I*N_CAC+cac]) + (_r[cac][0]*X[ORDER_S*N_CAC+cac]))*dt + - sqrt((_r[cac][3]*X[ORDER_I*N_CAC+cac]))*dem_sto__2[cac]- sqrt((_r[cac][2]*X[ORDER_I*N_CAC+cac]))*dem_sto__4[cac]+ sqrt((_r[cac][0]*X[ORDER_S*N_CAC+cac]))*dem_sto__1[cac];
        }
    }


    //TODO: drift of the diffusion
    //for(i=N_PAR_SV*N_CAC; i<(N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot); i++){
    //    f[i] = 0.0;
    //}

    /*automaticaly generated code:*/
    /*compute incidence:integral between t and t+1*/
    offset = N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot;

    
    o = 0;

    for (ts=0; ts<obs2ts[o]->n_ts_unique; ts++) {
        sum_inc = 0.0;
        for (n_cac=0; n_cac<obs2ts[o]->n_cac[ts]; n_cac++) {
            c = obs2ts[o]->cac[ts][n_cac][0];
            ac = obs2ts[o]->cac[ts][n_cac][1];
            cac = c*N_AC+ac;

            sum_inc += ( + (_r[cac][0]*X[ORDER_S*N_CAC+cac]))*dt + sqrt((_r[cac][0]*X[ORDER_S*N_CAC+cac]))*dem_sto__1[cac];
        }

        f[offset] = X[offset] +  sum_inc;
        offset++;
    }
    

    
    //y_pred (f) -> X (and we ensure that X is > 0.0)
    for(i=0; i<N_PAR_SV*N_CAC; i++){
        X[i] =  (f[i] < 0.0) ? 0.0 : f[i];
    }

    for(i=N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot; i<N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot +N_TS_INC_UNIQUE; i++){
        X[i] = (f[i] < 0.0) ? 0.0 : f[i];
    }
    

}


void step_sde_full(struct s_X *p_X, double t, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc)

{

    
    double *X = p_X->proj;
    double dt = p_X->dt;
    double *f = p_calc->y_pred;
    

    struct s_obs2ts **obs2ts = p_data->obs2ts;
    struct s_router **routers = p_data->routers;

    int i, c, ac, cac, n_cac, ts, o;
    double sum_inc = 0.0;
    int offset;
    
    double **par = p_par->natural;

    double _r[N_CAC][4];

    

    
    double dem_sto__0[N_CAC];
    double dem_sto__1[N_CAC];
    double dem_sto__2[N_CAC];
    double dem_sto__3[N_CAC];
    double dem_sto__4[N_CAC];
    double white_noise__0[N_CAC];

    


    for(cac=0;cac<N_CAC;cac++){
        


        

        
        _r[cac][0] = par[ORDER_r0][routers[ORDER_r0]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(X[ORDER_I*N_CAC+cac]+par[ORDER_i][routers[ORDER_i]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][1] = gsl_spline_eval(p_calc->spline[ORDER_mu_b][cac],t,p_calc->acc[ORDER_mu_b][cac])*gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][2] = gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac]);
        _r[cac][3] = (par[ORDER_v][routers[ORDER_v]->map[cac]]);

        
        dem_sto__0[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__1[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__2[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__3[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        dem_sto__4[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
        white_noise__0[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
    }

    for(c=0;c<N_C;c++) {
        for(ac=0; ac<N_AC; ac++) {
            cac = c*N_AC+ac;

            /*automaticaly generated code:*/
            /*ODE system*/
            
            f[0*N_CAC+cac] = X[0*N_CAC+cac] +  ( - (_r[cac][0]*X[ORDER_S*N_CAC+cac]) - (_r[cac][2]*X[ORDER_S*N_CAC+cac]) + (_r[cac][1]))*dt + - sqrt((_r[cac][0]*X[ORDER_S*N_CAC+cac]))*dem_sto__1[cac]- sqrt((_r[cac][2]*X[ORDER_S*N_CAC+cac]))*dem_sto__3[cac]+ sqrt((_r[cac][1]))*dem_sto__0[cac] + - (_r[cac][0]*X[ORDER_S*N_CAC+cac])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__0[cac];
            f[1*N_CAC+cac] = X[1*N_CAC+cac] +  ( - (_r[cac][3]*X[ORDER_I*N_CAC+cac]) - (_r[cac][2]*X[ORDER_I*N_CAC+cac]) + (_r[cac][0]*X[ORDER_S*N_CAC+cac]))*dt + - sqrt((_r[cac][3]*X[ORDER_I*N_CAC+cac]))*dem_sto__2[cac]- sqrt((_r[cac][2]*X[ORDER_I*N_CAC+cac]))*dem_sto__4[cac]+ sqrt((_r[cac][0]*X[ORDER_S*N_CAC+cac]))*dem_sto__1[cac] + + (_r[cac][0]*X[ORDER_S*N_CAC+cac])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__0[cac];
        }
    }


    //TODO: drift of the diffusion
    //for(i=N_PAR_SV*N_CAC; i<(N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot); i++){
    //    f[i] = 0.0;
    //}

    /*automaticaly generated code:*/
    /*compute incidence:integral between t and t+1*/
    offset = N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot;

    
    o = 0;

    for (ts=0; ts<obs2ts[o]->n_ts_unique; ts++) {
        sum_inc = 0.0;
        for (n_cac=0; n_cac<obs2ts[o]->n_cac[ts]; n_cac++) {
            c = obs2ts[o]->cac[ts][n_cac][0];
            ac = obs2ts[o]->cac[ts][n_cac][1];
            cac = c*N_AC+ac;

            sum_inc += ( + (_r[cac][0]*X[ORDER_S*N_CAC+cac]))*dt + + sqrt((_r[cac][0]*X[ORDER_S*N_CAC+cac]))*dem_sto__1[cac]  + + (_r[cac][0]*X[ORDER_S*N_CAC+cac])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__0[cac];
        }

        f[offset] = X[offset] +  sum_inc;
        offset++;
    }
    

    
    //y_pred (f) -> X (and we ensure that X is > 0.0)
    for(i=0; i<N_PAR_SV*N_CAC; i++){
        X[i] =  (f[i] < 0.0) ? 0.0 : f[i];
    }

    for(i=N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot; i<N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot +N_TS_INC_UNIQUE; i++){
        X[i] = (f[i] < 0.0) ? 0.0 : f[i];
    }
    

}


void step_sde_no_dem_sto(struct s_X *p_X, double t, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc)

{

    
    double *X = p_X->proj;
    double dt = p_X->dt;
    double *f = p_calc->y_pred;
    

    struct s_obs2ts **obs2ts = p_data->obs2ts;
    struct s_router **routers = p_data->routers;

    int i, c, ac, cac, n_cac, ts, o;
    double sum_inc = 0.0;
    int offset;
    
    double **par = p_par->natural;

    double _r[N_CAC][4];

    

    
    double white_noise__0[N_CAC];

    


    for(cac=0;cac<N_CAC;cac++){
        


        

        
        _r[cac][0] = par[ORDER_r0][routers[ORDER_r0]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(X[ORDER_I*N_CAC+cac]+par[ORDER_i][routers[ORDER_i]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][1] = gsl_spline_eval(p_calc->spline[ORDER_mu_b][cac],t,p_calc->acc[ORDER_mu_b][cac])*gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][2] = gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac]);
        _r[cac][3] = (par[ORDER_v][routers[ORDER_v]->map[cac]]);

        
        white_noise__0[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);
    }

    for(c=0;c<N_C;c++) {
        for(ac=0; ac<N_AC; ac++) {
            cac = c*N_AC+ac;

            /*automaticaly generated code:*/
            /*ODE system*/
            
            f[0*N_CAC+cac] = X[0*N_CAC+cac] +  ( - (_r[cac][0]*X[ORDER_S*N_CAC+cac]) - (_r[cac][2]*X[ORDER_S*N_CAC+cac]) + (_r[cac][1]))*dt + - (_r[cac][0]*X[ORDER_S*N_CAC+cac])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__0[cac];
            f[1*N_CAC+cac] = X[1*N_CAC+cac] +  ( - (_r[cac][3]*X[ORDER_I*N_CAC+cac]) - (_r[cac][2]*X[ORDER_I*N_CAC+cac]) + (_r[cac][0]*X[ORDER_S*N_CAC+cac]))*dt + + (_r[cac][0]*X[ORDER_S*N_CAC+cac])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__0[cac];
        }
    }


    //TODO: drift of the diffusion
    //for(i=N_PAR_SV*N_CAC; i<(N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot); i++){
    //    f[i] = 0.0;
    //}

    /*automaticaly generated code:*/
    /*compute incidence:integral between t and t+1*/
    offset = N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot;

    
    o = 0;

    for (ts=0; ts<obs2ts[o]->n_ts_unique; ts++) {
        sum_inc = 0.0;
        for (n_cac=0; n_cac<obs2ts[o]->n_cac[ts]; n_cac++) {
            c = obs2ts[o]->cac[ts][n_cac][0];
            ac = obs2ts[o]->cac[ts][n_cac][1];
            cac = c*N_AC+ac;

            sum_inc += ( + (_r[cac][0]*X[ORDER_S*N_CAC+cac]))*dt  + + (_r[cac][0]*X[ORDER_S*N_CAC+cac])*par[ORDER_sto][routers[ORDER_sto]->map[cac]]*white_noise__0[cac];
        }

        f[offset] = X[offset] +  sum_inc;
        offset++;
    }
    

    
    //y_pred (f) -> X (and we ensure that X is > 0.0)
    for(i=0; i<N_PAR_SV*N_CAC; i++){
        X[i] =  (f[i] < 0.0) ? 0.0 : f[i];
    }

    for(i=N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot; i<N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot +N_TS_INC_UNIQUE; i++){
        X[i] = (f[i] < 0.0) ? 0.0 : f[i];
    }
    

}


void step_sde_no_dem_sto_no_env_sto(struct s_X *p_X, double t, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc)

{

    
    double *X = p_X->proj;
    double dt = p_X->dt;
    double *f = p_calc->y_pred;
    

    struct s_obs2ts **obs2ts = p_data->obs2ts;
    struct s_router **routers = p_data->routers;

    int i, c, ac, cac, n_cac, ts, o;
    double sum_inc = 0.0;
    int offset;
    
    double **par = p_par->natural;

    double _r[N_CAC][4];

    

    

    


    for(cac=0;cac<N_CAC;cac++){
        


        

        
        _r[cac][0] = par[ORDER_r0][routers[ORDER_r0]->map[cac]]*par[ORDER_v][routers[ORDER_v]->map[cac]]*(X[ORDER_I*N_CAC+cac]+par[ORDER_i][routers[ORDER_i]->map[cac]])/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][1] = gsl_spline_eval(p_calc->spline[ORDER_mu_b][cac],t,p_calc->acc[ORDER_mu_b][cac])*gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][2] = gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac]);
        _r[cac][3] = (par[ORDER_v][routers[ORDER_v]->map[cac]]);

        
    }

    for(c=0;c<N_C;c++) {
        for(ac=0; ac<N_AC; ac++) {
            cac = c*N_AC+ac;

            /*automaticaly generated code:*/
            /*ODE system*/
            
            f[0*N_CAC+cac] = X[0*N_CAC+cac] +  ( - (_r[cac][0]*X[ORDER_S*N_CAC+cac]) - (_r[cac][2]*X[ORDER_S*N_CAC+cac]) + (_r[cac][1]))*dt;
            f[1*N_CAC+cac] = X[1*N_CAC+cac] +  ( - (_r[cac][3]*X[ORDER_I*N_CAC+cac]) - (_r[cac][2]*X[ORDER_I*N_CAC+cac]) + (_r[cac][0]*X[ORDER_S*N_CAC+cac]))*dt;
        }
    }


    //TODO: drift of the diffusion
    //for(i=N_PAR_SV*N_CAC; i<(N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot); i++){
    //    f[i] = 0.0;
    //}

    /*automaticaly generated code:*/
    /*compute incidence:integral between t and t+1*/
    offset = N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot;

    
    o = 0;

    for (ts=0; ts<obs2ts[o]->n_ts_unique; ts++) {
        sum_inc = 0.0;
        for (n_cac=0; n_cac<obs2ts[o]->n_cac[ts]; n_cac++) {
            c = obs2ts[o]->cac[ts][n_cac][0];
            ac = obs2ts[o]->cac[ts][n_cac][1];
            cac = c*N_AC+ac;

            sum_inc += ( + (_r[cac][0]*X[ORDER_S*N_CAC+cac]))*dt;
        }

        f[offset] = X[offset] +  sum_inc;
        offset++;
    }
    

    
    //y_pred (f) -> X (and we ensure that X is > 0.0)
    for(i=0; i<N_PAR_SV*N_CAC; i++){
        X[i] =  (f[i] < 0.0) ? 0.0 : f[i];
    }

    for(i=N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot; i<N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot +N_TS_INC_UNIQUE; i++){
        X[i] = (f[i] < 0.0) ? 0.0 : f[i];
    }
    

}



double likelihood(double x, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, const int ts, const int n, const double t)
{
    /*x is the predicted value from the model that we contrast with a time serie ts.
      Note: user should not use this function but get_log_likelihood
    */

    struct s_router **routers = p_data->routers;

    double like; /* likelihood value */

    double y = p_data->data[n][ts];

    double **par = p_par->natural;

    /*automaticaly generated code*/
    double gsl_mu = par[ORDER_rep][routers[ORDER_rep]->map[ts]]*x*gsl_spline_eval(p_calc->spline[ORDER_prop][ts],t,p_calc->acc[ORDER_prop][ts]);
    double gsl_sd = sqrt( pow(par[ORDER_phi][routers[ORDER_phi]->map[ts]],2)*pow(par[ORDER_rep][routers[ORDER_rep]->map[ts]],2)*pow(x,2)*pow(gsl_spline_eval(p_calc->spline[ORDER_prop][ts],t,p_calc->acc[ORDER_prop][ts]),2)+par[ORDER_rep][routers[ORDER_rep]->map[ts]]*x*gsl_spline_eval(p_calc->spline[ORDER_prop][ts],t,p_calc->acc[ORDER_prop][ts])*(-par[ORDER_rep][routers[ORDER_rep]->map[ts]]+1.0) );

    if (y > 0.0) {
        like=gsl_cdf_gaussian_P(y+0.5-gsl_mu, gsl_sd)-gsl_cdf_gaussian_P(y-0.5-gsl_mu, gsl_sd);
    } else {
        like=gsl_cdf_gaussian_P(y+0.5-gsl_mu, gsl_sd);
    }

    return sanitize_likelihood(like);
}

double obs_mean(double x, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, const int ts, const int n, const double t)
{
  /*x is the predicted value from the model that we contrast with a time serie ts*/
  struct s_router **routers = p_data->routers;

  double **par = p_par->natural;
  
  /*automaticaly generated code*/
  double mu = par[ORDER_rep][routers[ORDER_rep]->map[ts]]*x*gsl_spline_eval(p_calc->spline[ORDER_prop][ts],t,p_calc->acc[ORDER_prop][ts]);

  return mu;
}

double obs_var(double x, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, const int ts, const int n, const double t)
{
  /*x is the predicted value from the model that we contrast with a time serie ts*/
  struct s_router **routers = p_data->routers;

  double **par = p_par->natural;

  /*automaticaly generated code*/
  double var = pow(par[ORDER_phi][routers[ORDER_phi]->map[ts]],2)*pow(par[ORDER_rep][routers[ORDER_rep]->map[ts]],2)*pow(x,2)*pow(gsl_spline_eval(p_calc->spline[ORDER_prop][ts],t,p_calc->acc[ORDER_prop][ts]),2)+par[ORDER_rep][routers[ORDER_rep]->map[ts]]*x*gsl_spline_eval(p_calc->spline[ORDER_prop][ts],t,p_calc->acc[ORDER_prop][ts])*(-par[ORDER_rep][routers[ORDER_rep]->map[ts]]+1.0);

  return var;
}


double observation(double x, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, const int ts, const int n, const double t)
{
  /*x is the predicted value from the model that we contrast with a time serie ts*/
  struct s_router **routers = p_data->routers;

  double **par = p_par->natural;  

  /*return an observation of the process model*/

  /*automaticaly generated code*/
  double gsl_mu= par[ORDER_rep][routers[ORDER_rep]->map[ts]]*x*gsl_spline_eval(p_calc->spline[ORDER_prop][ts],t,p_calc->acc[ORDER_prop][ts]);
  double gsl_sd=sqrt(pow(par[ORDER_phi][routers[ORDER_phi]->map[ts]],2)*pow(par[ORDER_rep][routers[ORDER_rep]->map[ts]],2)*pow(x,2)*pow(gsl_spline_eval(p_calc->spline[ORDER_prop][ts],t,p_calc->acc[ORDER_prop][ts]),2)+par[ORDER_rep][routers[ORDER_rep]->map[ts]]*x*gsl_spline_eval(p_calc->spline[ORDER_prop][ts],t,p_calc->acc[ORDER_prop][ts])*(-par[ORDER_rep][routers[ORDER_rep]->map[ts]]+1.0));

  double yobs= gsl_mu+gsl_ran_gaussian(p_calc->randgsl, gsl_sd);

  return (yobs >0) ? yobs : 0.0;
}
