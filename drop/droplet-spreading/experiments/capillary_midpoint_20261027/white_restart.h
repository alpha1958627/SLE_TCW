/* Exact white-noise continuation only: the RNG is counter based.  Scratch
 * Poisson fields are rebuilt; height and cumulative diagnostics are retained.
 * Each continuation writes a NEW segment directory, never appends old logs. */
/* unistd declarations already enter through Basilisk grid/fpe.h; repeating
 * that system include here breaks this qcc parser's generated C. */
static int restart_offset=0, restart_total_steps=0, restart_stride=0;
static double white_lower=1e-6;
static const char *restart_path=NULL;
static double restart_now(void) { return (restart_offset+iter+1)*tcw_dt; }

typedef struct {
    char magic[16];
    unsigned long long seed;
    int n, completed, ranks, bridge;
    double dt, phi, lower, domain, floor;
    double initial_volume, initial_excess, minimum, previous_minimum;
    double floor_volume, pre_floor_minimum, joint_minimum;
    unsigned long long floor_cells;
    long solver_calls, solver_failures;
    int max_picard, max_cycles;
    double max_update, max_residual;
    double raw2[3], limited2[3], active[3], hits[3], removed2[3], minweight[3];
    int first[3];
} WhiteCheckpoint;

static void restart_read_header(WhiteCheckpoint *s) {
    FILE *fp=fopen(restart_path,"rb");
    if (!fp || fread(s,sizeof(*s),1,fp)!=1 || fgetc(fp)!=EOF) {
        fprintf(stderr,"CHECKPOINT_BAD_HEADER %s\n",restart_path); exit(41);
    }
    fclose(fp);
    if (memcmp(s->magic,"WHITE_CKPT_V1",13) || s->seed!=tcw_noise_seed ||
        s->n!=(1<<LEVEL_MIN) || s->ranks!=npe() || s->bridge!=tcw_bridge_depth ||
        s->dt!=tcw_dt || s->phi!=COEFF_VARPHI*tcw_varphi_scale ||
        s->lower!=white_lower || s->domain!=L0 || s->floor!=tcw_height_floor ||
        s->completed<=0 || s->completed>=restart_total_steps || s->solver_failures) {
        fprintf(stderr,"CHECKPOINT_CONFIG_MISMATCH\n"); exit(42);
    }
}

static void restart_load(void) {
    WhiteCheckpoint s; restart_read_header(&s);
    restart_offset=s.completed;
    char path[4096];
    int len=strlen(restart_path);
    if (len<6 || strcmp(restart_path+len-6,".state") || len>4000) exit(43);
    snprintf(path,sizeof(path),"%.*s.hfield",len-6,restart_path);
    int n=1<<LEVEL_MIN; size_t count=(size_t)n*n;
    double *values=malloc(count*sizeof(double)); int ok=1;
    if (pid()==0) {
        FILE *fp=fopen(path,"rb"); int nn=0;
        ok=fp && fread(&nn,sizeof(int),1,fp)==1 && nn==n &&
           fread(values,sizeof(double),count,fp)==count && fgetc(fp)==EOF;
        if (fp) fclose(fp);
    }
#if _MPI
    MPI_Bcast(&ok,1,MPI_INT,0,MPI_COMM_WORLD);
    if (ok) MPI_Bcast(values,n*n,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif
    if (!ok) { fprintf(stderr,"CHECKPOINT_BAD_FIELD\n"); exit(44); }
    foreach() {
        int ix=(int)floor((x-X0)/Delta), iy=(int)floor((y-Y0)/Delta);
        h[]=values[(size_t)iy*n+ix]; hi[]=1./h[];
    }
    free(values); boundary({h,hi});
    hmin=s.minimum; hmin1=s.previous_minimum;
    jcp_initial_volume=s.initial_volume; campaign_initial_excess=s.initial_excess;
    campaign_floor_volume=s.floor_volume; campaign_floor_cells=s.floor_cells;
    campaign_pre_floor_minimum=s.pre_floor_minimum; campaign_joint_raw_minimum=s.joint_minimum;
    tcw_solver_calls=s.solver_calls; tcw_solver_failures=s.solver_failures;
    tcw_solver_max_picard=s.max_picard; tcw_solver_max_mg_cycles=s.max_cycles;
    tcw_solver_max_update=s.max_update; tcw_solver_max_residual=s.max_residual;
    for (int k=0;k<3;k++) {
        campaign_flux_raw2[k]=s.raw2[k]; campaign_flux_limited2[k]=s.limited2[k];
        campaign_flux_active[k]=s.active[k]; campaign_flux_hits[k]=s.hits[k];
        scan_removed2[k]=s.removed2[k]; scan_minweight[k]=s.minweight[k]; scan_first[k]=s.first[k];
    }
    if (pid()==0) fprintf(stderr,"WHITE_RESTART completed=%d physical_time=%.17g next_noise_counter=%d floor_added=%.17g\n",restart_offset,restart_offset*tcw_dt,restart_offset,campaign_floor_volume);
}

static void restart_save(int completed) {
    WhiteCheckpoint s={0}; memcpy(s.magic,"WHITE_CKPT_V1",13);
    s.seed=tcw_noise_seed; s.n=1<<LEVEL_MIN; s.completed=completed;
    s.ranks=npe(); s.bridge=tcw_bridge_depth; s.dt=tcw_dt;
    s.phi=COEFF_VARPHI*tcw_varphi_scale; s.lower=white_lower;
    s.domain=L0; s.floor=tcw_height_floor;
    s.initial_volume=jcp_initial_volume; s.initial_excess=campaign_initial_excess;
    s.minimum=hmin; s.previous_minimum=hmin1;
    s.floor_volume=campaign_floor_volume; s.floor_cells=campaign_floor_cells;
    s.pre_floor_minimum=campaign_pre_floor_minimum; s.joint_minimum=campaign_joint_raw_minimum;
    s.solver_calls=tcw_solver_calls; s.solver_failures=tcw_solver_failures;
    s.max_picard=tcw_solver_max_picard; s.max_cycles=tcw_solver_max_mg_cycles;
    s.max_update=tcw_solver_max_update; s.max_residual=tcw_solver_max_residual;
    for (int k=0;k<3;k++) {
        s.raw2[k]=campaign_flux_raw2[k]; s.limited2[k]=campaign_flux_limited2[k];
        s.active[k]=campaign_flux_active[k]; s.hits[k]=campaign_flux_hits[k];
        s.removed2[k]=scan_removed2[k]; s.minweight[k]=scan_minweight[k]; s.first[k]=scan_first[k];
    }
    char stem[128], fieldname[160], temp[160], state[160];
    snprintf(stem,sizeof(stem),"checkpoints/step%09d",completed);
    snprintf(fieldname,sizeof(fieldname),"%s.hfield",stem);
    snprintf(temp,sizeof(temp),"%s.hfield.tmp",stem);
    snprintf(state,sizeof(state),"%s.state",stem);
    if (access(state,F_OK)==0 || access(fieldname,F_OK)==0) exit(45);
    write_hfield_binary(temp);
    if (pid()==0) {
        if (rename(temp,fieldname)) exit(46);
        snprintf(temp,sizeof(temp),"%s.state.tmp",stem);
        FILE *fp=fopen(temp,"wb");
        if (!fp || fwrite(&s,sizeof(s),1,fp)!=1 || fflush(fp) || fsync(fileno(fp))) exit(47);
        if (fclose(fp) || rename(temp,state)) exit(48);
        fprintf(stderr,"WHITE_CHECKPOINT completed=%d physical_time=%.17g state=%s field=%s\n",completed,completed*tcw_dt,state,fieldname);
    }
#if _MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}
