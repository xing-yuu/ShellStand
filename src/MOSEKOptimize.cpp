#include "MOSEKOptimize.h"
std::vector < cinolib::vec3d>* msk::get_poly_center(cinolib::DrawableTrimesh<> m) {
    std::vector<std::vector<uint>> faces = m.vector_polys();
    std::vector<cinolib::vec3d> vertexs = m.vector_verts();
    std::vector < cinolib::vec3d>* poly_center = new std::vector < cinolib::vec3d>;
    for (uint i = 0; i < faces.size(); i++) {
        cinolib::vec3d vi = vertexs[faces[i][0]];
        cinolib::vec3d ui = vertexs[faces[i][1]];
        cinolib::vec3d wi = vertexs[faces[i][2]];
        poly_center->push_back((vi + ui + wi) / 3);
    }
    return poly_center;
}
std::vector < std::vector < uint>>* msk::poly_link_poly(const cinolib::DrawableTrimesh<> m) {
    std::vector < std::vector < uint>>* poly_link = new std::vector < std::vector < uint>>;
    std::vector<std::vector<uint>> faces = m.vector_polys();
    for (uint pid = 0; pid < faces.size(); pid++) {
        std::set<uint> s;
        std::vector < uint> l;
        for (uint vid = 0; vid < 3; vid++) {
            std::vector<uint> pl = m.vert_ordered_polys_star(faces[pid][vid]);
            for (uint k = 0; k < pl.size(); k++)
                if (pl[k] < pid)
                    s.insert(pl[k]);
        }
        for (std::set<uint>::iterator it = s.begin(); it != s.end(); it++)
            l.push_back(*it);
        poly_link->push_back(l);
    }
    return poly_link;
}



std::vector<double>* msk::MOSEK_deltam(std::vector<double> m0, std::vector<double> rho, cinolib::vec3d target_center, std::vector<double> upper_bound, double lower_bound, cinolib::DrawableTrimesh<> m) {
    std::vector<std::vector<uint>> faces = m.vector_polys();
    std::vector<cinolib::vec3d> vertexs = m.vector_verts();
    std::vector<cinolib::vec3d>* poly_center = msk::get_poly_center(m);
    const double rho_coef = 1, smooth_coef = 1;
    const double CONERROR = 0.01;
    

    double cx = 0, cy = 0;
    for (uint i = 0; i < faces.size(); i++) {
        cx += ((poly_center->at(i)[0] - target_center[0]) * m0[i] * m.poly_area(i));
        cy += ((poly_center->at(i)[1] - target_center[1]) * m0[i] * m.poly_area(i));
    }
    MSKtask_t          task = NULL;
    MSKrescodee        r = MSK_RES_OK;
    MSKint32t          i, j;
    MSKenv_t           env = NULL;
    MSKint32t          numvar = faces.size(),
                       numcon = 2;



    r = MSK_makeenv(&env, NULL);
    if (r == MSK_RES_OK) {
        /* Create the optimization task. */
        r = MSK_maketask(NULL, numcon, numvar, &task);
        if (r == MSK_RES_OK)
            r = MSK_linkfunctotaskstream(task, MSK_STREAM_LOG, NULL, printstr);
        /* Append 'numcon' empty constraints. The constraints will initially have no bounds. */
        if (r == MSK_RES_OK)
            r = MSK_appendcons(task, numcon);

        /* Append 'numvar' variables. The variables will initially be fixed at zero (x=0). */
        if (r == MSK_RES_OK)
            r = MSK_appendvars(task, numvar);

        if (r == MSK_RES_OK)
            r = MSK_putcfix(task, 0.0);


        //Constraint condition
        MSKboundkeye bkc[] = { MSK_BK_RA,MSK_BK_RA };
        double blc[] = { -cx - fabs(cx * CONERROR), -cy - fabs(cy * CONERROR) };
        double buc[] = { -cx + fabs(cx * CONERROR), -cy + fabs(cy * CONERROR) };

        MSKboundkeye* bkx = new MSKboundkeye[int(numvar)];
        double* blx = new double[int(numvar)];
        double* bux = new double[int(numvar)];
        for (j = 0; j < numvar && r == MSK_RES_OK; ++j) {
            bkx[j]=(MSK_BK_RA);
            blx[j]=lower_bound - m0[int(j)];
            bux[j]= std::max(upper_bound[int(j)], lower_bound) - m0[int(j)];
        }
        for (j = 0; j < numvar && r == MSK_RES_OK; ++j) {
            /* Set the bounds on variable j. blx[j] <= x_j <= bux[j] */
            if (r == MSK_RES_OK)
                r = MSK_putcj(task, j, 0);

            if (r == MSK_RES_OK)
                r = MSK_putvarbound(task,
                    j,           /* Index of variable.*/
                    bkx[j],      /* Bound key.*/
                    blx[j],      /* Numerical value of lower bound.*/
                    bux[j]);     /* Numerical value of upper bound.*/
               
                
            // Input  A
            if (r == MSK_RES_OK)
                r = MSK_putaij(task,
                    0,
                    j,
                    (poly_center->at(j)[0] - target_center[0]) * m.poly_area(j));
            if (r == MSK_RES_OK)
                r = MSK_putaij(task,
                    1,
                    j,
                    (poly_center->at(j)[1] - target_center[1]) * m.poly_area(j));
                    
        }
        /* Set the bounds on constraints.
       for i=1, ...,numcon : blc[i] <= constraint i <= buc[i] */
        for (i = 0; i < numcon && r == MSK_RES_OK; ++i)

            r = MSK_putconbound(task,
                i,           /* Index of constraint.*/
                bkc[i],      /* Bound key.*/
                blc[i],      /* Numerical value of lower bound.*/
                buc[i]);     /* Numerical value of upper bound.*/


        //OBJ condition/smooth
        std::vector < std::vector < uint>>* poly_link = msk::poly_link_poly(m);
        std::vector < uint> numvaild(numvar);
        for (j = 0; j < numvar && r == MSK_RES_OK; ++j) {
            for (uint pid = 0; pid < poly_link->at(uint(j)).size() && r == MSK_RES_OK; pid++) {
                numvaild[uint(j)]++;
                numvaild[poly_link->at(uint(j))[pid]]++;
                r = MSK_putqobjij(
                            task,
                            MSKint32t(j),
                            MSKint32t(poly_link->at(uint(j))[pid]),
                            -smooth_coef);
            }
        }
        for (j = 0; j < numvar && r == MSK_RES_OK; ++j) {
            r = MSK_putqobjij(
                task,
                j,
                j,
                rho[uint(j)] * rho_coef + numvaild[uint(j)] * smooth_coef);
        }


        /* Minimize objective function. */
        if (r == MSK_RES_OK)
            r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE);

 
        if (r == MSK_RES_OK) {

            MSKrescodee trmcode;

            /* Run optimizer */
            r = MSK_optimizetrm(task, &trmcode);

            /* Print a summary containing information
               about the solution for debugging purposes. */
            MSK_solutionsummary(task, MSK_STREAM_LOG);

            if (r == MSK_RES_OK)
            {
                MSKsolstae solsta;

                if (r == MSK_RES_OK)
                    r = MSK_getsolsta(task,
                        MSK_SOL_ITR,
                        &solsta);
                switch (solsta)
                {
                case MSK_SOL_STA_OPTIMAL:
                {
                    double *xx = (double*)calloc(numvar, sizeof(double));
                    if (xx)
                    {
                        MSK_getxx(task,
                            MSK_SOL_ITR,    /* Request the basic solution. */
                            xx);
                        std::vector<double>* res = new std::vector<double>(xx, xx + numvar);
                        return res;
                    }
                    else {
                        r = MSK_RES_ERR_SPACE;
                        std::string vi = "ERROR.obj";
                        m.save(vi.data());
                    }
                    break;
                }
                case MSK_SOL_STA_DUAL_INFEAS_CER:
                case MSK_SOL_STA_PRIM_INFEAS_CER:
                    printf("Primal or dual infeasibility certificate found.\n");
                    break;
                case MSK_SOL_STA_UNKNOWN:
                {
                    char symname[MSK_MAX_STR_LEN];
                    char desc[MSK_MAX_STR_LEN];

                    /* If the solutions status is unknown, print the termination code
                       indicating why the optimizer terminated prematurely. */

                    MSK_getcodedesc(trmcode,
                        symname,
                        desc);

                    printf("The solution status is unknown.\n");
                    printf("The optimizer terminitated with code: %s\n", symname);
                    break;
                }
                default:
                    printf("Other solution status.\n");
                    break;
                }
            }
        }

        if (r != MSK_RES_OK)
        {
            /* In case of an error print error code and description. */
            char symname[MSK_MAX_STR_LEN];
            char desc[MSK_MAX_STR_LEN];

            printf("An error occurred while optimizing.\n");
            MSK_getcodedesc(r,
                symname,
                desc);
            printf("Error %s - '%s'\n", symname, desc);
        }

        /* Delete the task and the associated data. */
        MSK_deletetask(&task);
    }
    MSK_deleteenv(&env);
}














void msk::Mosek_test() {
    uint NUMCON = 2;
    uint NUMVAR = 30000;

    double        c[] = { 0.0, -1.0, 0.0 };

    MSKboundkeye  bkc[] = { MSK_BK_LO ,MSK_BK_LO };
    double        blc[] = { 1.0,1.0 };
    double        buc[] = { +MSK_INFINITY ,+MSK_INFINITY };

    MSKboundkeye  bkx[] = { MSK_BK_LO,
                           MSK_BK_LO,
                           MSK_BK_LO
    };
    double        blx[] = { 0.0,
                           0.0,
                           0.0
    };
    double        bux[] = { +MSK_INFINITY,
                            +MSK_INFINITY,
                            +MSK_INFINITY
    };

    MSKint32t     aptrb[] = { 0,   1,   2 },
        aptre[] = { 1,   2,   3 },
        asub[] = { 0,   0,   0 };
    double        aval[] = { 1.0, 1.0, 1.0 };




    MSKint32t     i, j;
    double        *xx=new double [NUMVAR];

    MSKenv_t      env = NULL;
    MSKtask_t     task = NULL;
    MSKrescodee   r;

    /* Create the mosek environment. */
    r = MSK_makeenv(&env, NULL);

    if (r == MSK_RES_OK)
    {
        /* Create the optimization task. */
        r = MSK_maketask(env, NUMCON, NUMVAR, &task);

        if (r == MSK_RES_OK)
        {
            r = MSK_linkfunctotaskstream(task, MSK_STREAM_LOG, NULL, printstr);

            /* Append 'NUMCON' empty constraints.
             The constraints will initially have no bounds. */
            if (r == MSK_RES_OK)
                r = MSK_appendcons(task, NUMCON);

            /* Append 'NUMVAR' variables.
             The variables will initially be fixed at zero (x=0). */
            if (r == MSK_RES_OK)
                r = MSK_appendvars(task, NUMVAR);

            /* Optionally add a constant term to the objective. */
            if (r == MSK_RES_OK)
                r = MSK_putcfix(task, 0.0);
            for (j = 0; j < NUMVAR && r == MSK_RES_OK; ++j)
            {
                /* Set the linear term c_j in the objective.*/
                if (r == MSK_RES_OK)
                    r = MSK_putcj(task, j, 0);

                /* Set the bounds on variable j.
                 blx[j] <= x_j <= bux[j] */
                if (r == MSK_RES_OK)
                    r = MSK_putvarbound(task,
                        j,           /* Index of variable.*/
                        MSK_BK_LO,      /* Bound key.*/
                        0,      /* Numerical value of lower bound.*/
                        MSK_INFINITY);     /* Numerical value of upper bound.*/

  // Input column j of A */
                if (r == MSK_RES_OK)
                    r = MSK_putaij(task, 0, j, 1);
                if (r == MSK_RES_OK)
                    r = MSK_putaij(task, 1, j, 1);
                //r = MSK_putacol(task,

                  //  j,                 /* Variable (column) index.*/
             //       aptre[j] - aptrb[j], /* Number of non-zeros in column j.*/
            //        asub + aptrb[j],   /* Pointer to row indexes of column j.*/
             //       aval + aptrb[j]);  /* Pointer to Values of column j.*/

            }

            /* Set the bounds on constraints.
               for i=1, ...,NUMCON : blc[i] <= constraint i <= buc[i] */
            for (i = 0; i < NUMCON && r == MSK_RES_OK; ++i)
                r = MSK_putconbound(task,
                    i,           /* Index of constraint.*/
                    bkc[i],      /* Bound key.*/
                    blc[i],      /* Numerical value of lower bound.*/
                    buc[i]);     /* Numerical value of upper bound.*/
            for (j = 0; j < NUMVAR && r == MSK_RES_OK; ++j) {
                if (r == MSK_RES_OK)
                    r = MSK_putqobjij(
                        task,
                        MSKint32t(j),
                        MSKint32t(j),
                        1);
            }
            
            if (r == MSK_RES_OK)
            {
                MSKrescodee trmcode;

                /* Run optimizer */
                r = MSK_optimizetrm(task, &trmcode);

                /* Print a summary containing information
                   about the solution for debugging purposes*/
                MSK_solutionsummary(task, MSK_STREAM_MSG);

                if (r == MSK_RES_OK)
                {
                    MSKsolstae solsta;
                    int j;

                    MSK_getsolsta(task, MSK_SOL_ITR, &solsta);

                    switch (solsta)
                    {
                    case MSK_SOL_STA_OPTIMAL:
                        MSK_getxx(task,
                            MSK_SOL_ITR,    /* Request the interior solution. */
                            xx);

                        printf("Optimal primal solution\n");
                        //      for (j = 0; j < NUMVAR; ++j)
                         //         printf("x[%d]: %e\n", j, xx[j]);

                        break;

                    case MSK_SOL_STA_DUAL_INFEAS_CER:
                    case MSK_SOL_STA_PRIM_INFEAS_CER:
                        printf("Primal or dual infeasibility certificate found.\n");
                        break;

                    case MSK_SOL_STA_UNKNOWN:
                        printf("The status of the solution could not be determined. Termination code: %d.\n", trmcode);
                        break;

                    default:
                        printf("Other solution status.");
                        break;
                    }
                }
                else
                {
                    printf("Error while optimizing.\n");
                }
            }

            if (r != MSK_RES_OK)
            {
                /* In case of an error print error code and description. */
                char symname[MSK_MAX_STR_LEN];
                char desc[MSK_MAX_STR_LEN];

                printf("An error occurred while optimizing.\n");
                MSK_getcodedesc(r,
                    symname,
                    desc);
                printf("Error %s - '%s'\n", symname, desc);
            }
        }
        MSK_deletetask(&task);
    }
    MSK_deleteenv(&env);

}

