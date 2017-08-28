/*
* This file constitutes a part of RNANR program.
*
* RNANR is free software to which the terms of the GNU Affero General 
* Public License apply as published by the Free Software Foundation, 
* either of version 3, or (at your option) any later version.
*
* This program is distributed without ANY WARRANTY or RESPONSIBILITY.
* You can distribute and change its contents, youŕe not allowed to SELL it.
* For more precisions, see the GNU Affero General Public License.
*
* You should have received a copy of the GNU Affero General Public License
* along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

/**************************************************/
/*                   counting                    */
/************************************************/

/********************************************/
/*    operators for the locopt structures  */
/******************************************/

// computes Boltzmann factor and passes its value by reference to first parameter
void apply_energy_term(TYPE * returned_value, double dG, TYPE val){ 
  // modifying term for energy (Boltzmann factor)
  //printf("reformatos: %f -> %f\n", dG, exp(-dG/RT));
  TYPE dG_cast_neg, RT_cast, Boltzmann_factor;
  INIT(dG_cast_neg);
  INIT(RT_cast);
  INIT(Boltzmann_factor);
  SET_TYPE_VAL_FROM_DOUBLE(dG_cast_neg, -dG); // casts -dG to TYPE
  SET_TYPE_VAL_FROM_DOUBLE(RT_cast, RT); // casts RT to TYPE
  
  /* This is unfortunately necessary because mpfr does only have assignment operations.
   * The original expression is 'exp(-dG/RT)*val' : [starts here] =>
   */
  SET_TYPE_VAL(Boltzmann_factor, dG_cast_neg);
  DIVIDE(Boltzmann_factor, RT_cast);
  EXPONENT(Boltzmann_factor);
  MULTIPLY(Boltzmann_factor, val);
  // <= [ends here] 
  //printf("reformatos: %f -> %f\n", dG, exp(-dG/RT));
  SET_TYPE_VAL(*returned_value, Boltzmann_factor);   //Boltzmann's cte already contains temperature
  // cleaning 
  CLEAR(dG_cast_neg); 
  CLEAR(RT_cast); 
  CLEAR(Boltzmann_factor);
}

/***********   end utilities     ************/

/* Creates a table of Boltzmann partition*/
void init_partition_table(int n){
  int x, y; 
  partition_function_table=(TYPE **) malloc ((n+1)*sizeof(TYPE *));
  for (x=0; x<=n; x++){
    partition_function_table[x]=(TYPE *) malloc ((n+1)*sizeof(TYPE));
    for (y=0; y<=n; y++){
	  INIT(partition_function_table[x][y]);	
      SET_TYPE_VAL_FROM_DOUBLE(partition_function_table[x][y], 0.);/* yann */
    }
  } 
}


// changes pointers within table (modifies directly flat_table via pointers)
void rewire(int x, int y, linktab * linktable){
	if(linktable->length > 0){
		flat_table[x][y] = linktable->links[0].fcell;
		for(int	i = 0; i<(linktable->length-1); i++){
			linktable->links[i].fcell->next = linktable->links[i+1].fcell;
		}
		linktable->links[linktable->length-1].fcell->next = NULL;
	}
}


void get_partition_function(TYPE * partition_function, plain_sequence * rna){
  TYPE partition_fci; /* terms of partition function */
  TYPE intermediate; /* stores intermediate values for partition function */
  INIT(partition_fci);
  INIT(intermediate);
  
  int x, y, i, j, i0; 
  flat_cell * current_flat_structure;
  int current_base_pair;
  int thickness; 
  int alpha; /* number of helixes in loop */
  int beta; /* number of unpaired bases */
  char *rna_seq = (char*) malloc(sizeof(char)*(rna->size+1));
  strncpy(rna_seq, rna->label+1, rna->size);
  rna_seq[rna->size] = '\0';
  
  init_partition_table(rna->size);
  vrna_md_t md;
  set_model_details(&md);
  md.temperature = TEMP;  
  E_fold_cp = vrna_fold_compound(rna_seq, &md, VRNA_OPTION_DEFAULT);
  
  for (x=rna->size; x>=1; x--){
    for (y=x+1; y<=rna->size; y++){
	   linktab * linktable = start_linktab(); // create a table of links to be sorted after
      if (flat_table[x][y]!=NULL){
	current_flat_structure=flat_table[x][y];
	do {
	  current_base_pair=current_flat_structure->current;
	  i=get_flat_structure_start(current_base_pair);
	  j=get_flat_structure_end(current_base_pair); 
      if ((i==x) && (j==y) && !is_entirely_contained(i,j,rna)){ /* exception to cases where i=1 OR j=rna->size */
	    /* thickness= 1 */
	    apply_energy_term(&partition_fci, stacking_energy(x-1,y+1), partition_function_table[x+1][y-1]); /* yann */
        add_link_element(linktable, DOUBLE_CAST(partition_fci), current_flat_structure);
	  }
	  else{
		SET_TYPE_VAL_FROM_DOUBLE(partition_fci, 1.);  
        alpha = 0;
        beta = 0;
        i0 = x - 1;
		int *ptypes =  (int *) malloc (sizeof(int)); /* table of pair types for multibranch loops */
		int *si = (int *) malloc (sizeof(int));
		int *sj = (int *) malloc (sizeof(int));
		ptypes[0] =  get_type(rna -> label[y+1], rna -> label[x-1]); /* closing pair type */ 
		si[0] = y;
		sj[0] = x;     
	    while(current_base_pair != 0){
		  alpha ++;
	      i=get_flat_structure_start(current_base_pair);
	      j=get_flat_structure_end(current_base_pair);
	      /* external loop contribution */
		  if(is_entirely_contained(x,y,rna)){
			 apply_energy_term(&partition_fci, ext_loop_energy(i, j, rna), partition_fci); 
		  }
	      get_flat_structure_suffix(current_base_pair);
	      ptypes =  (int *) realloc (ptypes, (alpha+1)*sizeof(int));
	      si = (int *) realloc (si, (alpha+1)*sizeof(int));
	      sj = (int *) realloc (sj, (alpha+1)*sizeof(int));	      
	      ptypes[alpha] = get_type(rna -> label[i], rna -> label[j]); /* enclosed pair type(s) */
	      si[alpha] = i-1;
	      sj[alpha] = j+1;
		  beta += (i - i0 - 1);  
	      i0 = j;
	      thickness=get_BP(i,j);
	      if (thickness>MIN_HELIX_LENGTH) 
		thickness=MIN_HELIX_LENGTH;
		  apply_energy_term(&intermediate, helix_energy(i,j,thickness), partition_function_table[i+thickness][j-thickness]);
		  MULTIPLY(partition_fci, intermediate);
	      current_base_pair = get_flat_structure_suffix(current_base_pair);
	    } /* end while */
	    if((x != 1) || (y != rna->size)){  
	      if(alpha == 1){ /* internal loop or bulge since stack was taken care of before*/
			apply_energy_term(&partition_fci, internal_loop_energy(x-1, i, j, y+1, rna), partition_fci);  
		  }else{ /*multiloop*/
			beta += (y - i0);
			apply_energy_term(&partition_fci, multiloop_energy(alpha, beta, ptypes, si, sj), partition_fci);
		  }	
		}
		free(ptypes);
		free(si);
		free(sj);
	  }/* endif */
	  /* next line is sum because we are adding energies for different secondary structures*/
	  ADD(partition_function_table[x][y], partition_fci);
	  add_link_element(linktable, DOUBLE_CAST(partition_fci), current_flat_structure); // adds link to the existing link table
	  current_flat_structure=current_flat_structure->next;
	}while (current_flat_structure != NULL);
      }/*endif else */
      else if ((x>1) && (y<rna->size)){ /* exception to cases where i=1 OR j=rna->size */
		SET_TYPE_VAL_FROM_DOUBLE(intermediate, 1.);  
		apply_energy_term(&partition_function_table[x][y], hairpin_energy(x-1,y+1, rna), intermediate);
	  }
	  linktable = sort(linktable);
	  rewire(x, y, linktable);
	  free_linktab(linktable);
    }/* end for y */
  }/* end for x */
   
  SET_TYPE_VAL(* partition_function, partition_function_table[1][rna->size]);
  
  // freeing memory
  CLEAR(partition_fci);
  CLEAR(intermediate);
  free(rna_seq); 
}

/* probabilistic backtracking of random sampling (redundant and non-redundant) */


void stochastic_backtrack_locally_optimal_structure_rec(int x, int y, plain_sequence * rna, char * structure, 
 TYPE * energy, tr_node ** actual_node, TYPE * denomin);

/* DP_t and struc_count are references to execution time of DP and number of structures respectively */
folding* stochastic_backtrack_locally_optimal_structures(int number_of_structures, plain_sequence * rna, 
 int is_non_redun, int use_timer, float zqlimit, int * struc_count, double * DP_t){
  /*clock_t begin_time = clock();
  clock_t current_time;
  double time_spent;*/
  /* Precompute/cache #locOpts */
  
  int i,j, continuing_on;
  long int max_struct;
  TYPE total_Boltzmann;
  TYPE cumulative_parf;
  TYPE coverage; // % of Boltzmann factor that is covered by currently returned solution
  TYPE* denomin  = (TYPE*) malloc (sizeof(TYPE));  // denominator used to ponderate the weight of forbidden structures
  INIT(total_Boltzmann);
  INIT(cumulative_parf);
  SET_TYPE_VAL_FROM_DOUBLE(cumulative_parf, 0.);
  INIT(coverage);
  INIT(*denomin);
  
  // intiates RNG
  INIT_RNG();
  
  /* timer-related variables */
  struct timespec start_time, end_time;
  
  /* struct holding RNA structure and energy*/
  folding *folding_energy = (folding*)malloc (sizeof(folding));
  max_struct = count_all_locally_optimal_structures(rna);
  if((max_struct < number_of_structures) && !(is_non_redun)){
	 printf("Error: maximum number of structures is %ld. Use number lower or equal to this with option -s.\n\
	 Alternatively, you can also use option -z for non-redundant sampling. \n\n", max_struct); 
	 exit(EXIT_SUCCESS);
  } 
  
  /* creates tree to ensure non-redundant sampling */
  tr_node * non_red_tree = create_root();
  tr_node * actual_node = non_red_tree; // pointer to actual node
    
  get_partition_function(&total_Boltzmann, rna); // total Boltzmann will be set to value of partition function computed by this fci
  folding_energy->structures = (char**) malloc (sizeof(char*));
  folding_energy->part_fcis = (double*) malloc (sizeof(double));
  folding_energy->energy_ref = (float*) malloc (sizeof(float));
  if(use_timer) clock_gettime(CLOCK_MONOTONIC, &start_time);
  /* this is basically fancy 'for' loop, since it allows two different types of samplings without much redundancy */
  i = 1;
  continuing_on = 1;
  while(continuing_on){
	//printf("\nGenerating new structure: %d\n", i);
	folding_energy->structures = (char**) realloc (folding_energy->structures, (i+1)*sizeof(char*));
    folding_energy->part_fcis = (double*) realloc (folding_energy->part_fcis, (i+1)*sizeof(TYPE));
    folding_energy->energy_ref = (float*) realloc (folding_energy->energy_ref, (i+1)*sizeof(float));
	folding_energy->part_fcis[i] = 1.;
	
    char* structure;
    TYPE* partit_fci  = (TYPE*) malloc (sizeof(TYPE));
    INIT(*partit_fci);
    structure = (char*) malloc ((rna->size+1)*sizeof(char));
    for (j=1;j<=rna->size;j++){
      structure[j] = '.';
    }
    SET_TYPE_VAL_FROM_DOUBLE(*partit_fci, 1.);
    SET_TYPE_VAL(*denomin, partition_function_table[1][rna->size]); // = Z
    stochastic_backtrack_locally_optimal_structure_rec(1, rna->size, rna, structure, partit_fci, &actual_node, denomin);   
    folding_energy->structures[i] = structure;
    folding_energy->part_fcis[i] = DOUBLE_CAST(*partit_fci);
    ADD(cumulative_parf, *partit_fci);
    //printf("Total sum of chosen partition part: %0.3f \n", cumulative_parf);
    folding_energy->energy_ref[i] = get_reference_energy(E_fold_cp, structure, rna);
    if(!is_non_redun) actual_node = traceback_to_root(actual_node, *partit_fci);
    /*current_time = clock();
    time_spent = (double)(current_time-begin_time) / CLOCKS_PER_SEC;
    printf("clock: - %f - ", time_spent);
    print_structure(structure, *energy, get_reference_energy(E_fold_cp, structure, rna), rna);*/
    CLEAR(*partit_fci);
    free(partit_fci);
    i+=1;
    if(zqlimit != 0){
		SET_TYPE_VAL(coverage, cumulative_parf);
		DIVIDE(coverage, total_Boltzmann);
		if((double)zqlimit <= DOUBLE_CAST(coverage)){
			continuing_on = 0;
			*struc_count = i-1;
		}
	}else{
		if(i > number_of_structures){
			continuing_on = 0;
		}	
	}
  } 
  if(use_timer){
    clock_gettime(CLOCK_MONOTONIC, &end_time);
    *DP_t = (double)(end_time.tv_sec - start_time.tv_sec) + (double)(end_time.tv_nsec - start_time.tv_nsec)/1.0e9 ;
  } 
  CLEAR(total_Boltzmann);
  CLEAR(cumulative_parf);
  CLEAR(coverage);
  CLEAR(*denomin);
  return folding_energy;
}

void stochastic_backtrack_flat_structure_rec(flat_cell * current_flat_structure, plain_sequence * rna, 
 char * structure, TYPE * partit_fci, tr_node ** actual_node, TYPE * denomin){
  int current_base_pair;
  int i, j; 
  int thickness; 
						
  current_base_pair = current_flat_structure->current;
  while(current_base_pair != 0){
    i=get_flat_structure_start(current_base_pair);
    j=get_flat_structure_end(current_base_pair);
	thickness=get_BP(i,j); 
    if (thickness>MIN_HELIX_LENGTH) 
      thickness=MIN_HELIX_LENGTH;
    add_helix(i, j, thickness, structure);
    stochastic_backtrack_locally_optimal_structure_rec(i+thickness, j-thickness, rna, structure, partit_fci, actual_node, denomin);
    current_base_pair=get_flat_structure_suffix(current_base_pair);
  } /* end while */  
}

void stochastic_backtrack_locally_optimal_structure_rec(int x, int y, plain_sequence * rna, char * structure, 
 TYPE * partit_fci, tr_node ** actual_node, TYPE * denomin){
  double e_contribution; /* e_contribution - contribution of rna element (ML etc.) to energy */ 	 
  // init of all TYPE
  TYPE r, zero;   // zero is necessary to be used in comparison
  TYPE partition_fci, local_part_weight; // temporary vars to store intermediate computations in
  TYPE intermediate1, intermediate2, weighted_prob; // intermediates to compute Boltzmann factor of ALL AVAILABLE structures; first two also used later
  /* intermediate1 and intermediate2 are used multiple times to speed up algorithm. However, every time they are used, they're set
   * just before (using fcis SET_TYPE_VAL, SET_TYPE_VAL_FROM_DOUBLE or apply_energy_term)
   */
  
  INIT(r);
  INIT(zero);
  SET_TYPE_VAL_FROM_DOUBLE(zero, 0.);
  INIT(partition_fci);
  INIT(local_part_weight);
  INIT(intermediate1);
  INIT(intermediate2);
  INIT(weighted_prob);
  
  //initializaton of all other types
  int i, j, i0;
  flat_cell * current_flat_structure;
  int current_base_pair;
  int thickness; 
  int alpha; /*number of helixes in pair*/
  int beta; /*number of unpaired bases*/
  
  SET_TYPE_VAL(intermediate1, partition_function_table[x][y]);
  DIVIDE(intermediate1, *denomin);
  total_weight_par(&intermediate2, *actual_node);
  MULTIPLY(intermediate2, intermediate1);
  
  //intermediate1 = partition_function_table[x][y]/(*denomin);
  //intermediate2 = total_weight_par(*actual_node)*intermediate1;
  SET_TYPE_VAL(weighted_prob, partition_function_table[x][y]);
  SUBSTRACT(weighted_prob, intermediate2); 
  RANDOM_NUMBER(r, weighted_prob);
  
  // clear intermediate values which aren't needed anymore
  CLEAR(weighted_prob);
  
  if (flat_table[x][y]!=NULL){
  	current_flat_structure=flat_table[x][y];
  	do {
  	  current_base_pair=current_flat_structure->current;
  	  i=get_flat_structure_start(current_base_pair);
  	  j=get_flat_structure_end(current_base_pair); 
      if ((i==x) && (j==y) && !( (i==1) && (j==rna->size) )){
		  e_contribution = stacking_energy(x-1,y+1);
		  /* r = r - (apply_energy_term(e_contribution, partition_function_table[x+1][y-1] - 
		   *  (tr_node_weight(current_flat_structure->current, *actual_node) * partition_function_table[x][y])/(*denomin))
		   * <=>
		   * r = r + (tr_node_weight(current_flat_structure->current, *actual_node) * partition_function_table[x][y])/(*denomin)) - 
		   * (apply_energy_term(e_contribution, partition_function_table[x+1][y-1]
		   * 
		   *  ^ This is interpretation used in following lines
		   * starts here => */ 
		  tr_node_weight(&intermediate1, current_flat_structure->current, *actual_node); // intermediate1 = weigth of forbidden structures on actual node in this context
		  // intemediate1 gets reused to reduce execution time, ditto to 2
		  MULTIPLY(intermediate1, partition_function_table[x][y]);     
		  DIVIDE(intermediate1, *denomin);								
		  apply_energy_term(&intermediate2, e_contribution, partition_function_table[x+1][y-1]);
		  SUBSTRACT(intermediate1, intermediate2);
		  ADD(r, intermediate1);
		  // <= ends here
          if (IS_LOWER(r,zero)){
			add_base_pair(x,y,structure);
			*actual_node = add_if_nexists(current_flat_structure->current, *actual_node);
			apply_energy_term(partit_fci, e_contribution, *partit_fci);
			// reusing intermediate1 and intermediate2 in another context to speed up algorithm
			/* *denomin *= exp(e_contribution/RT)*partition_function_table[x+1][y-1]/partition_function_table[x][y]) 
			 * (starts here) => */
			SET_TYPE_VAL(intermediate2, partition_function_table[x+1][y-1]);
			DIVIDE(intermediate2, partition_function_table[x][y]); 
			apply_energy_term(&intermediate1, e_contribution, intermediate2);
			MULTIPLY(*denomin, intermediate1);
			// <= ends here
			stochastic_backtrack_locally_optimal_structure_rec(x+1, y-1, rna, structure, partit_fci, actual_node, denomin);
			
			CLEAR(r);
			CLEAR(zero);
			CLEAR(partition_fci);
			CLEAR(local_part_weight);
			CLEAR(intermediate1);
			CLEAR(intermediate2);
			return ;
          }
  	  }
  	  else{
		SET_TYPE_VAL_FROM_DOUBLE(local_part_weight, 1.);
		SET_TYPE_VAL_FROM_DOUBLE(partition_fci, 1.);
		alpha = 0;
		beta = 0;
		i0 = x - 1;
		int *ptypes = (int *) malloc (sizeof(int)); /* table of pair types for multibranch loops */
		int *si = (int *) malloc (sizeof(int));
		int *sj = (int *) malloc (sizeof(int));
		ptypes[0] =  get_type(rna -> label[y+1], rna -> label[x-1]); /* closing pair type */ 
		si[0] = y;
		sj[0] = x;
  	    while(current_base_pair != 0){
		  alpha ++;
  	      i=get_flat_structure_start(current_base_pair);
  	      j=get_flat_structure_end(current_base_pair);
  	      ptypes =  (int *) realloc (ptypes, (alpha+1)*sizeof(int));
  	      si = (int *) realloc (si, (alpha+1)*sizeof(int));
	      sj = (int *) realloc (sj, (alpha+1)*sizeof(int));	      
	      ptypes[alpha] = get_type(rna -> label[i], rna -> label[j]); /* enclosed pair type(s) */
	      si[alpha] = i-1;
	      sj[alpha] = j+1;
		  beta += (i - i0 - 1);  
          i0 = j;
  	      thickness=get_BP(i,j); 
  	      if (thickness>MIN_HELIX_LENGTH) 
  		      thickness=MIN_HELIX_LENGTH;
  		  e_contribution = helix_energy(i,j,thickness);
  		  apply_energy_term(&intermediate1, e_contribution, partition_function_table[i+thickness][j-thickness]);
  		  MULTIPLY(partition_fci, intermediate1); /* yann */
  	      apply_energy_term(&local_part_weight, e_contribution, local_part_weight);
  	      // apply bonuses/penalties for external loop
  	      if(is_entirely_contained(x, y, rna)){
			apply_energy_term(&partition_fci, ext_loop_energy(i, j, rna), partition_fci);  
			apply_energy_term(&local_part_weight, ext_loop_energy(i, j, rna), local_part_weight);
          }
  	      current_base_pair=get_flat_structure_suffix(current_base_pair);
  	    } /* end while */
  	    if((x != 1) || (y != rna->size)){  
	      if(alpha == 1){ /* internal loop or bulge since stack was taken care of before*/		  
			e_contribution =  internal_loop_energy(x-1, i, j, y+1, rna);
		  }else{ /*multiloop*/
			beta += (y - i0);
			e_contribution =  multiloop_energy(alpha, beta, ptypes, si, sj);
		  }
		  apply_energy_term(&partition_fci, e_contribution, partition_fci);	
          apply_energy_term(&local_part_weight, e_contribution, local_part_weight);
		}
  	    /* Now we have the overall "weight" of flat structure, should we choose it? 
  	     * r = r - (partition_fci - (tr_node_weight(&intermediate1, current_flat_structure->current, *actual_node) * partition_function_table[x][y]/ *denomin)
  	     *  <=>
  	     * r = r + (tr_node_weight(&intermediate1, current_flat_structure->current, *actual_node) * partition_function_table[x][y]/ *denomin - partition_fci)
  	     * ^ used interpretation
  	     * (starts here) => */
  	    tr_node_weight(&intermediate1, current_flat_structure->current, *actual_node);
  	    MULTIPLY(intermediate1, partition_function_table[x][y]);
		DIVIDE(intermediate1, *denomin);
		SUBSTRACT(intermediate1, partition_fci);
		ADD(r, intermediate1);
  	    // <= (ends here)
        if (IS_LOWER(r, zero)){
			 MULTIPLY(*partit_fci, local_part_weight);
			 // reusing intermediate1 in different context
			 /* *denomin *= (partition_fci/partition_function_table[x][y]);
			  * (starts here) =>*/
			 SET_TYPE_VAL(intermediate1, partition_fci);
			 DIVIDE(intermediate1, partition_function_table[x][y]);
			 MULTIPLY(*denomin, intermediate1);
			 // <= (ends here)
			 *actual_node = add_if_nexists(current_flat_structure->current, *actual_node);
             stochastic_backtrack_flat_structure_rec(current_flat_structure, rna, structure, partit_fci, actual_node, denomin);
             
             CLEAR(r);
             CLEAR(zero);
             CLEAR(partition_fci);
			 CLEAR(local_part_weight);
			 CLEAR(intermediate1);
			 CLEAR(intermediate2);
			 /* Warning: Potential memory leak below */
             return;
          }
          free(ptypes);
          free(si);
          free(sj); 	    
  	  }/* endif */
  	  current_flat_structure=current_flat_structure->next;
  	} while (current_flat_structure != NULL);
  }/*endif else */
  else{
	if((x != 1) || (y != rna->size))
		apply_energy_term(partit_fci, hairpin_energy(x-1, y+1, rna), *partit_fci);
		
	CLEAR(r);
	CLEAR(zero);
	CLEAR(partition_fci);
	CLEAR(local_part_weight);
	CLEAR(intermediate1);
	CLEAR(intermediate2);
    return;
  }
}
