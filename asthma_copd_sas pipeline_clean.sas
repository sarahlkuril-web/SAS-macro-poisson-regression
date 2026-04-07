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
        pov_pct        = INPUT(pov_pct,        best12.);
        public_ins_pct = INPUT(public_ins_pct,  best12.);
        PM25_Mean      = INPUT(PM25_Mean,        best12.);
    RUN;

%mend import_subset;


/* Import all 4 subsets */
%import_subset(sheet=Asthma_IP, out=AIP);
%import_subset(sheet=Asthma_OP, out=AOP);   /* FIX: removed space from sheet name */
%import_subset(sheet=COPD_IP,   out=CIP);   /* FIX: removed space from sheet name */
%import_subset(sheet=COPD_OP,   out=COP);   /* FIX: removed space from sheet name */


/*----- STEP 2: MACRO FOR POISSON GLM -----*/
/* FIX: Added ODS OUTPUT to capture parameter estimates for rr_summary */

%macro run_poisson(data=, outcome=, offset=, covars=, outname=);

    ODS OUTPUT ParameterEstimates = est_&outname;

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

%mend run_poisson;


/*----- STEP 3: MACRO FOR GLMM WITH RANDOM MONTH EFFECT -----*/
/* FIX: Parameterized so it loops across all four subsets;
        equivalent to glmmPQL in R (PROC GLIMMIX, method RSPL) */

%macro run_glmm(data=, outcome=, offset=, covars=, outname=);

    ODS OUTPUT ParameterEstimates = est_glmm_&outname;

    PROC GLIMMIX DATA=&data METHOD=RSPL;
        CLASS month;
        MODEL &outcome = &covars / DIST=POISSON LINK=LOG OFFSET=&offset SOLUTION;
        RANDOM INTERCEPT / SUBJECT=month;
        OUTPUT OUT=glmm_&outname RESID=residuals PRED=predicted;
    RUN;

    ODS OUTPUT CLOSE;

%mend run_glmm;


/*----- STEP 4: RUN GLM AND GLMM ACROSS ALL FOUR SUBSETS -----*/

%macro run_all_subsets;

    %let subsets = AIP AOP CIP COP;

    %do i = 1 %to 4;
        %let ds = %scan(&subsets, &i);

        /* --- Poisson GLM --- */
        %run_poisson(
            data    = &ds,
            outcome = cases,
            offset  = log_pop,
            covars  = over64_pct pov_pct college_pct public_ins_pct tempmax_avg PM25_90pct,
            outname = glm_&ds
        );

        /* --- GLMM with random month effect --- */
        %run_glmm(
            data    = &ds,
            outcome = cases,
            offset  = log_pop,
            covars  = over64_pct pov_pct college_pct public_ins_pct tempmax_avg PM25_90pct,
            outname = &ds
        );

    %end;

%mend run_all_subsets;

%run_all_subsets;


/*----- STEP 5: MORAN'S I (spatial autocorrelation check) -----*/
/* Calls PROC IML to compute Moran's I from residuals + weights matrix.
   FIX: Corrected S0 formula (sum of all elements of W);
        added loop so this can be called per-subset residual dataset. */

%macro morans_i(resid_data=, resid_var=, weights_csv=);

    /* Load spatial weights matrix */
    PROC IMPORT DATAFILE="&weights_csv"
        OUT=wmat
        DBMS=CSV REPLACE;
    RUN;

    PROC IML;

        /* Load residuals */
        USE &resid_data;
            READ ALL VAR {&resid_var} INTO e;
        CLOSE &resid_data;

        /* Load W matrix (square; rows = cols = n ZCTAs) */
        USE wmat;
            READ ALL INTO W;
        CLOSE wmat;

        n  = NROW(e);
        S0 = SUM(W);           /* FIX: SUM(W) gives total of all weights */

        ebar = e[:];
        z    = e - ebar;

        num = z` * W * z;      /* FIX: added transpose (z`) for conformability */
        den = z` * z;

        I = (n / S0) * (num / den);
        PRINT "Moran's I statistic:" I;

    QUIT;

%mend morans_i;

/* Example calls — one per subset residual dataset */
%morans_i(
    resid_data  = glm_AIP,
    resid_var   = residuals,
    weights_csv = C:\Users\vivre\Documents\Dissertation\wmat.csv
);

%morans_i(
    resid_data  = glm_AOP,
    resid_var   = residuals,
    weights_csv = C:\Users\vivre\Documents\Dissertation\wmat.csv
);

%morans_i(
    resid_data  = glm_CIP,
    resid_var   = residuals,
    weights_csv = C:\Users\vivre\Documents\Dissertation\wmat.csv
);

%morans_i(
    resid_data  = glm_COP,
    resid_var   = residuals,
    weights_csv = C:\Users\vivre\Documents\Dissertation\wmat.csv
);


/*----- STEP 6: SUMMARIZE MODEL OUTPUT (Relative Risk + 95% CI) -----*/
/* Works for parameter estimate datasets captured via ODS OUTPUT
   from either PROC GENMOD or PROC GLIMMIX.
   FIX: References est_glm_<subset> datasets created in Step 2,
        not a hardcoded "glm_estimates" dataset. */

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
        param_data = est_glm_&ds,
        beta_var   = Estimate,
        multiplier = 1,
        outname    = rr_raw_&ds
    );

    /* RR per 10-unit increase (e.g. PM2.5, temperature) */
    %rr_summary(
        param_data = est_glm_&ds,
        beta_var   = Estimate,
        multiplier = 10,
        outname    = rr_x10_&ds
    );

%end;
