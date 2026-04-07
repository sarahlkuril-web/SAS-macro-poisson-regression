/*----------------------------------------------------
File Name: Asthma/COPD Exacerbation Analysis Pipeline
Author:    Sarah L. Kuril
Project:   Dissertation
Description: Analytic macros for Poisson Regression GLM
and GLMM across four encounter subsets:
    - Asthma Inpatient  (AIP)
    - Asthma Outpatient (AOP)
    - COPD Inpatient    (CIP)
    - COPD Outpatient   (COP)
------------------------------------------------------*/


/*----- STEP 1: LOAD DATA -----*/

%macro import_subset(sheet=, out=);

    PROC IMPORT DATAFILE="C:\Users\vivre\Documents\Dissertation\pop_subset_rearranged.xlsx"
        OUT=&out
        DBMS=XLSX REPLACE;
        SHEET="&sheet";
    RUN;

    /* Ensure correct variable formats */
    DATA &out;
        SET &out;
        zcta5          = PUT(zcta5, 5.);
        month          = PUT(month, 2.);
        pov_pct        = INPUT(pov_pct, best12.);
        public_ins_pct = INPUT(public_ins_pct, best12.);
        PM25_Mean      = INPUT(PM25_Mean, best12.);
    RUN;

%mend import_subset;


/* Import all 4 subsets */
%import_subset(sheet=Asthma_IP, out=AIP);
%import_subset(sheet=Asthma_OP, out=AOP);
%import_subset(sheet=COPD_IP,   out=CIP); 
%import_subset(sheet=COPD_OP,   out=COP);


/*----- STEP 2: MACRO FOR POISSON GLM -----*/

%macro run_poisson_glm(data=, outcome=, offset=, covars=, outname=);

    ODS OUTPUT ParameterEstimates = est_&outname; /* ODS OUTPUT to capture parameter estimates for rr_summary */

    PROC GENMOD DATA=&data;
        MODEL &outcome = &covars / DIST=POISSON LINK=LOG OFFSET=&offset;
        OUTPUT OUT=&outname RESRAW=residuals PRED=predicted;
    RUN;

    ODS OUTPUT CLOSE;

    TITLE "Poisson Regression Residuals: &outname";
    PROC PRINT DATA=&outname (OBS=20);
        VAR zcta5 &outcome predicted residuals;
    RUN;
    TITLE;   /* Clear title after each print */

%mend run_poisson_glm;


/*----- STEP 3: MACRO FOR GLMM WITH RANDOM MONTH EFFECT -----*/
/* Macro loops across all four subsets */

%macro run_poisson_glmm(data=, outcome=, offset=, covars=, outname=);

    ODS OUTPUT ParameterEstimates = est_glmm_&outname;

    PROC GLIMMIX DATA=&data METHOD=RSPL; /* RSPL=Residual Pseudo-Likelihood, a default estimation method for non-normal data in GLMM */
        CLASS month;
        MODEL &outcome = &covars / DIST=POISSON LINK=LOG OFFSET=&offset SOLUTION;
        RANDOM INTERCEPT / SUBJECT=month;
        OUTPUT OUT=glmm_&outname RESID=residuals PRED=predicted;
    RUN;

    ODS OUTPUT CLOSE;

%mend run_poisson_glmm;


/*----- STEP 4: RUN GLM AND GLMM ACROSS ALL FOUR SUBSETS -----*/

%macro run_all_subsets;

    %let subsets = AIP AOP CIP COP;

    %do i = 1 %to 4;
        %let ds = %scan(&subsets, &i);

        /* --- Poisson GLM --- */
        %run_poisson_glm(
            data    = &ds,
            outcome = cases,
            offset  = log_pop,
            covars  = over64_pct pov_pct college_pct public_ins_pct tempmax_avg PM25_90pct,
            outname = glm_&ds
        );

        /* --- GLMM with random month effect --- */
        %run_poisson_glmm(
            data    = &ds,
            outcome = cases,
            offset  = log_pop,
            covars  = over64_pct pov_pct college_pct public_ins_pct tempmax_avg PM25_90pct,
            outname = &ds
        );

    %end;

%mend run_all_subsets;

%run_all_subsets;

/*-------------------------------------------------------------------------------------*/
/*---------RECC END HERE, FOR STATISTICAL DIVERGENCE REASONS AS NOTED BELOW -----------*/
/*-------------------------------------------------------------------------------------*/



/*----- STEP 5: MORAN'S I (spatial autocorrelation check) -----*/
/* Calls PROC IML to compute Moran's I from residuals + pre-exported weights matrix from a CSV (k=10, great-circle distance)*/

%macro morans_i(resid_data=, resid_var=, weights_csv=);

    /* Load spatial weights (W) matrix */
    /* This reads a CSV file containing pre-exported weights matrix using the 10 nearest neighbors (k=10, great-circle distance)*/
    PROC IMPORT DATAFILE="&weights_csv"
        OUT=wmat
        DBMS=CSV REPLACE;
    RUN;

    PROC IML;

        /* Load residuals */
        /* ZCTA-averaged residuals from the GLM */
        USE &resid_data;
            READ ALL VAR {&resid_var} INTO e;
        CLOSE &resid_data;

        /* Load W matrix (square; rows = cols = n ZCTAs) */
        USE wmat;
            READ ALL INTO W;
        CLOSE wmat;

        n  = NROW(e);
        S0 = SUM(W);           /* SUM(W) gives total of all weights */

        ebar = e[:];           /* Grand mean of residuals (scalar) */
        z    = e - ebar;

        num = z` * W * z;      /* Transpose (z`) added for conformability
                                  Numerator: spatial lag of mean-centered residuals. z` is the row transpose needed for matrix conformability (1 x n * n x n * n x 1 = scalar). */
        den = z` * z;          /* Denominator: total sum of squares of z. */

        I = (n / S0) * (num / den);  /* Moran's I statistic. */
        PRINT "Moran's I statistic:" I;

    QUIT;

%mend morans_i;

/* Run Moran's I for each subset using its GLM residual dataset. */

%morans_i(
    resid_data  = glm_AIP,
    resid_var   = residuals,
    weights_csv = C:\Users\vivre\Documents\Dissertation\wmat.csv
);
/* AIP = Asthma Inpatient. Residuals sourced from glm_AIP output
   dataset created in Step 2. Confirm row count matches n ZCTAs
   in wmat.csv before interpreting I. */

%morans_i(
    resid_data  = glm_AOP,
    resid_var   = residuals,
    weights_csv = C:\Users\vivre\Documents\Dissertation\wmat.csv
);
/* AOP = Asthma Outpatient. Same wmat assumed — revisit if the
   outpatient ZCTA coverage differs from inpatient. */

%morans_i(
    resid_data  = glm_CIP,
    resid_var   = residuals,
    weights_csv = C:\Users\vivre\Documents\Dissertation\wmat.csv
);
/* CIP = COPD Inpatient. */

%morans_i(
    resid_data  = glm_COP,
    resid_var   = residuals,
    weights_csv = C:\Users\vivre\Documents\Dissertation\wmat.csv
);
/* COP = COPD Outpatient. */


/*----- STEP 6: SUMMARIZE MODEL OUTPUT (Relative Risk + 95% CI) -----*/
/* Works for parameter estimate datasets captured via ODS OUTPUT from either PROC GENMOD or PROC GLIMMIX. */

%macro rr_summary(param_data=, beta_var=, multiplier=1, outname=);

    DATA &outname;
        SET &param_data;
        RR     = EXP(&beta_var * &multiplier);
        RR_LCL = EXP(LowerCL   * &multiplier);
        RR_UCL = EXP(UpperCL   * &multiplier);
        LABEL
            RR     = "Relative Risk"
            RR_LCL = "95% LCL"
            RR_UCL = "95% UCL";
    RUN;

    TITLE "Relative Risk Summary: &outname (multiplier=&multiplier)";
    PROC PRINT DATA=&outname LABEL NOOBS;
        VAR Parameter RR RR_LCL RR_UCL;
    RUN;
    TITLE;

%mend rr_summary;


/* --- RR summaries for GLM estimates across all four subsets --- */

%let subsets = AIP AOP CIP COP;

%do i = 1 %to 4;
    %let ds = %scan(&subsets, &i);

    /* Raw RR (per 1-unit increase) */
    %rr_summary(
        param_data = est_glm_&ds,    /* references est_glm_<subset> datasets created in Step 2 */
        beta_var   = Estimate,
        multiplier = 1,
        outname    = rr_raw_&ds
    );

    /* RR per 10-unit increase (e.g. PM2.5, temperature) */
    %rr_summary(
        param_data = est_glm_&ds,  /* references est_glm_<subset> datasets created in Step 2 */
        beta_var   = Estimate,
        multiplier = 10,
        outname    = rr_x10_&ds
    );

%end;

/*NOTES: The most important methodological note to keep in mind: the R code constructs the weights matrix 
on the fly using knearneigh(k=10, longlat=TRUE), which uses great-circle (Haversine) distances — appropriate 
for lat/lon coordinates. The SAS version simply reads whatever CSV you give it, so the validity of the 
comparison depends entirely on that CSV having been exported from the identical kNN specification in R via 
something like listw2mat(lstwaip). If there's any mismatch in k, distance metric, or ZCTA ordering between 
the two, the I statistics will diverge for reasons unrelated to the underlying spatial pattern. 
Therefore, will run GLMM + Bayesian in R instead for analysis*/
