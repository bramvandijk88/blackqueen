MersenneTwister = require('mersennetwister');
const yargs = require('yargs');
const fs = require('fs');

if (typeof window == "undefined") {
	Simulation = require('./cacatooM.js')
//	All = require('/Volumes/Beta_Structural/SBQ_Modelling/JavaScript/cacatoo/cacatoo-master/lib/all.js') 
}

const { hideBin } = require('yargs/helpers')
const argv = yargs(hideBin(process.argv)).argv
//console.log("argv", argv)


        let num_CG = 3				// Number of CGs in system
        let global_mut = 0.000005		//	Chance for each gene in each cell to mutate into non-production per timestep
        let global_death = 0.1		// Chance of death for each cell per timestep
        let base_repro_chance = 0.6		// Default chance a non-producer may be able to reproduce
        let starting_omni_proportion = 1.0
        let cost_per_CG = 0.05		// Proportionate penalty for each CG 
        let CG_radius = 1		// radius of Moore neighbourhood for benefits from common goods
        let temp_var
        var count_alive=0			// How many cells are alive
        let births = 0						// Count num births
        let kills = 0 
        let maxsteps = 50000                             // How many time steps the model continues to run            
        let num_col = 200                                   // Number of columns (width of your grid)
        let num_row = 200								// Number of columns (width of your grid)
        let diffuse = 1				// If >0, then margolis diffusion will be enabled
        let write_out_number = 1000		// How many entries are desired in the output files
        let write_out_interval = Math.round(maxsteps / write_out_number)	// convert he amount of writeouts to be proportional to duration of the simulation
        
		let Chromosome_HGT_rate = 0		// Chance of a CG gene in a living cell to copy itself over its ortholog in a neighbouring cell's genome
		let Zero_HGT_ratio = 1		// Bias in HGT of non-producing alleles. 1 = zeros transfer at same rate as 1 (producing gene); 0 = zeros never transfer (treats them as deleted over inactivated)
        
        let filn_suffix

        var list = [0,1]
        let permCounts = undefined
        let masterPerm = undefined
        let run_num = 1


if (typeof argv.num_CG !== "undefined") {
	num_CG = argv.num_CG
	console.log("num_CG -- user specified", num_CG)
	} else {
		console.log("num_CG\t", num_CG)
	}

if (typeof argv.global_mut !== "undefined") {
	global_mut = argv.global_mut
	console.log("global_mut -- user specified\t", global_mut)
} else {
	console.log("global_mut\t", global_mut)
}

if (typeof argv.global_death !== "undefined") {
	global_death = argv.global_death
	console.log("global_death -- user specified\t", global_death)
} else {
	console.log("global_death\t", global_death)
}

if (typeof argv.cost_per_CG !== "undefined") {
	cost_per_CG = argv.cost_per_CG
	console.log("cost_per_CG -- user specified\t", cost_per_CG)
} else {
	console.log("cost_per_CG\t", cost_per_CG)
}

if (typeof argv.base_repro_chance !== "undefined") {
	base_repro_chance = argv.base_repro_chance
	console.log("base_repro_chance -- user specified\t", base_repro_chance)
} else {
	console.log("base_repro_chance\t", base_repro_chance)
}


if (typeof argv.CG_radius !== "undefined") {	
	CG_radius = argv.CG_radius
	console.log("CG_radius -- user specified\t", CG_radius)
	} else {
	console.log("CG_radius\t", CG_radius)
}

if (typeof argv.maxsteps !== "undefined") {	
	maxsteps = argv.maxsteps
	console.log("maxsteps -- user specified\t", maxsteps)
	} else {
	console.log("maxsteps\t", maxsteps)
}

if (typeof argv.run_num !== "undefined") {	
	run_num = argv.run_num
	console.log("run_num -- user specified\t", run_num)
	} else {
	console.log("run_num\t", run_num)
}

if (typeof argv.diffuse !== "undefined") {	
	diffuse = argv.diffuse
	console.log("diffuse -- user specified\t", diffuse)
	} else {
	console.log("diffuse\t", diffuse)
}

if (typeof argv.write_out_interval !== "undefined") {	
	write_out_interval = argv.write_out_interval
	console.log("write_out_interval -- user specified\t", write_out_interval)
	} else {
	console.log("write_out_interval\t", write_out_interval)
}

if (typeof argv.write_out_number !== "undefined") {	
	write_out_number = argv.write_out_number
	console.log("write_out_number -- user specified\t", write_out_number)
	} else {
	console.log("write_out_number\t", write_out_number)
}

if (typeof argv.Chromosome_HGT_rate !== "undefined") {
	Chromosome_HGT_rate = argv.Chromosome_HGT_rate
	console.log("Chromosome_HGT_rate -- user specified", Chromosome_HGT_rate)
	} else {
		console.log("Chromosome_HGT_rate\t", Chromosome_HGT_rate)
	}
	
if (typeof argv.Zero_HGT_ratio !== "undefined") {
	Zero_HGT_ratio = argv.Zero_HGT_ratio
	console.log("Zero_HGT_ratio -- user specified", Zero_HGT_ratio)
	} else {
		console.log("Zero_HGT_ratio\t", Zero_HGT_ratio)
	}


console.log("")


var dir = './data_out_'+run_num
if (!fs.existsSync(dir)){
    fs.mkdirSync(dir);
} 

filn_suffix = CG_radius+"_"+global_mut+"_"+cost_per_CG+"_"+global_death+"_"+num_CG+"_"+base_repro_chance+"_"+Chromosome_HGT_rate+"_"+Zero_HGT_ratio

let filn = dir+"/Run_Info-"+filn_suffix+".note"
if(fs.existsSync(filn)){
	fs.unlinkSync(filn)
}
fs.appendFileSync(filn, "run_num\t"+run_num+"\n"+"num_CG\t"+num_CG+"\n"+"global_mut\t"+global_mut+"\n"+"global_death\t"+global_death+"\n")
fs.appendFileSync(filn, "base_repro_chance\t"+base_repro_chance+"\n"+"starting_omni_proportion\t"+starting_omni_proportion+"\n")
fs.appendFileSync(filn, "cost_per_CG\t"+cost_per_CG+"\n"+"CG_radius\t"+CG_radius+"\n"+"maxsteps\t"+maxsteps+"\n")
fs.appendFileSync(filn, "num_col\t"+num_col+"\n"+"num_row\t"+num_row+"\n"+"margolis diffusion\t"+diffuse+"\n")
fs.appendFileSync(filn, "write_out_interval\t"+write_out_interval+"\n")
fs.appendFileSync(filn, "Chromosome_HGT_rate\t"+Chromosome_HGT_rate+"\n")
fs.appendFileSync(filn, "Zero_HGT_ratio\t"+Zero_HGT_ratio+"\n")

filn = dir+"/count_alive_"+filn_suffix+".txt"
if(fs.existsSync(filn)){
	fs.unlinkSync(filn)
}

filn = dir+"/production_"+filn_suffix+".txt"	
if(fs.existsSync(filn)){
	fs.unlinkSync(filn)
}

filn = dir+"/kill_count_"+filn_suffix+".txt"	
if(fs.existsSync(filn)){
	fs.unlinkSync(filn)
}

filn = dir+"/CheaterCounts_"+filn_suffix+".txt"	
if(fs.existsSync(filn)){
	fs.unlinkSync(filn)
}

grid_dir = dir+"/GridStates_"+filn_suffix+"/"	
if (!fs.existsSync(grid_dir)){
    fs.mkdirSync(grid_dir);
} 

filn = dir+"/CheaterCounts_"+"doppleganger_"+filn_suffix+".txt"	
if(fs.existsSync(filn)){
	fs.unlinkSync(filn)
}

grid_dir = dir+"/GridStates_"+"doppleganger_"+filn_suffix+"/"	
if (!fs.existsSync(grid_dir)){
    fs.mkdirSync(grid_dir);
} 

filn = dir+"/HGT_debug_"+filn_suffix+".txt"	
if(fs.existsSync(filn)){
//	fs.unlinkSync(filn)
}

// Re-calc write out interval in case the user provided a new number of desired write outs
write_out_interval = Math.round(maxsteps / write_out_number)

    /*-----------------------Start user-defined code ---------------------*/

   let sim; // Declare a variable named "sim" globally, so that we can access our cacatoo-simulation from wherever we need. 

    /**
    * function cacatoo() contains all the user-defined parts of a cacatoo-model. Configuration, update rules, what is displayed or plotted, etc. It's all here.
    */
    function cacatoo() {
        /*
            1. SETUP. 
        */
        let config =
        {
            title: "sBQ",                        // The name of your cacatoo-simulation
            description: "",         // description, if you needed
            maxtime: maxsteps,                             // How many time steps the model continues to run            
            ncol: num_col,                                   // Number of columns (width of your grid)
            nrow: num_row,		                              // Number of rows (height of your grid)
            wrap: [true, true],                         // Wrapped boundary conditions? [COLS, ROWS]   
            scale: 2,				                      // Scale of the grid (nxn pixels per grid point)
        }

        /*
            1. SETUP. (continued) 
        */
        sim = new Simulation(config)                          // Initialise the Cacatoo simulation
        sim.makeGridmodel("cells")                              // Build a new Gridmodel within the simulation called "cells"
		sim.makeGridmodel("dopple")   						//	For invasions, or for doubling the runs per instance - create a parallel grid 

        /*
            2. DEFINING THE RULES. Below, the user defines the nextState function. This function will be applied for each grid point when we will update the grid later. 
        */

		//	Initialise a grid - write an initial accounting through reportGenomes() function
        sim.cells.initialise = function () {
            sim.initialGrid(sim.cells, "alive", 0, 1.0)

            for (let i = 0; i < sim.cells.nc; i++) for (let j = 0; j < sim.cells.nr; j++) {
                if (sim.cells.rng.random() < starting_omni_proportion) {
                    sim.cells.grid[i][j].alive = 1
                    sim.cells.grid[i][j].genome = new Array(num_CG).fill(1)             
                    temp_var = sim.cells.grid[i][j].genome
                    sim.cells.grid[i][j].production = temp_var.reduce((a, b) => a + b, 0)
                }
            }
            
            sim.cells.reportGenomes()
        }


		//	For invasions, or for doubling the runs per instance, create a parallel grid - write an initial accounting through reportGenomes() function
        sim.dopple.initialise = function () {
            sim.initialGrid(sim.dopple, "alive", 0, 1.0)

            for (let i = 0; i < sim.dopple.nc; i++) for (let j = 0; j < sim.dopple.nr; j++) {
                if (sim.dopple.rng.random() < starting_omni_proportion) {
                    sim.dopple.grid[i][j].alive = 1
                    sim.dopple.grid[i][j].genome = new Array(num_CG).fill(1)             
                    temp_var = sim.dopple.grid[i][j].genome
                    sim.dopple.grid[i][j].production = temp_var.reduce((a, b) => a + b, 0)
                }
            }
            sim.dopple.reportGenomes()
        }
 

		//	For doppleganger grid - count how many cells have each possible genotype        
        sim.dopple.reportGenomes = function() {
            makePermCounts()
            for (let col_iter = 0; col_iter < sim.dopple.nc; col_iter++) {
                for (let row_iter = 0; row_iter < sim.dopple.nr; row_iter++) {
                   if (this.grid[col_iter][row_iter].alive == 1) {
                        let t_genome = this.grid[col_iter][row_iter].genome.toString()
                        t_val = permCounts.get(t_genome)
                        permCounts.set(t_genome, t_val+1)
                    }
                }

            }
        }                


		// randomise an array - given the list of living neighbours and randomised so there is nothing deterministic about growth order.
        function shuffle(array) {
            for (let q = array.length - 1; q > 0; q--) {
                let r = Math.floor(Math.random() * (q + 1));
                [array[q], array[r]] = [array[r], array[q]];
            }
            return array
        }


       //Borrowed from the net. Create every combo of [0,1] for the number of CGs in the model 
       getPermutations = function(list, maxLen) {
        // Copy initial values as arrays
            var perm = list.map(function(val) {
                return [val];
            });
            // Our permutation generator
            var generate = function(perm, maxLen, currLen) {
            // Reached desired length
                if (currLen === maxLen) {
                    return perm;
                }
                // For each existing permutation
                for (var i = 0, len = perm.length; i < len; i++) {
                    var currPerm = perm.shift();
                    // Create new permutation
                    for (var k = 0; k < list.length; k++) {
                        perm.push(currPerm.concat(list[k]));
                    }
                }
                // Recurse
                return generate(perm, maxLen, currLen + 1);
            };
            // Start with size 1 because of initial values
            return generate(perm, maxLen, 1);
        }

        
        //	Create every genotype permutation from the number of CGs in sim
        makePermCounts = function() {
            permCounts = new Map();
            var res = getPermutations(list, num_CG)
            for (r = 0; r < res.length; r++) {
                t_string = res[r].toString()
                permCounts.set(t_string, 0)
            }
        }

        
        //	Initial creation of every genotype permutation from the number of CGs in sim 
        First_makePermCounts = function() {
            permCounts = new Map();
            var res = getPermutations(list, num_CG)
            for (r = 0; r < res.length; r++) {
                t_string = res[r].toString()
                permCounts.set(t_string, 0)
            }
            masterPerm = permCounts
        }
  
        
        //	Count how many cells have each possible genotype  
        sim.cells.reportGenomes = function() {
            makePermCounts()
            for (let col_iter = 0; col_iter < sim.cells.nc; col_iter++) {
                for (let row_iter = 0; row_iter < sim.cells.nr; row_iter++) {
                   if (this.grid[col_iter][row_iter].alive == 1) {
                        let t_genome = this.grid[col_iter][row_iter].genome.toString()
                        t_val = permCounts.get(t_genome)
                        permCounts.set(t_genome, t_val+1)
                    }
                }
            }
        }                
        
		
		//	feed a reference gridpoint. Return all living cells within a box within the radius of CG effect in this sim.
        sim.cells.findLivingNeighbours = function (i,j) {
            let living_neighbours = [] // empty array where you will store grid points in
            for(let ii= (i - CG_radius); ii<= (i + CG_radius); ii++)
            {
                for(let jj= (j - CG_radius); jj<= (j + CG_radius); jj++)
                {
                    if(sim.cells.getGridpoint(ii,jj) == undefined) {}
                    else if(sim.cells.getGridpoint(ii,jj).alive == 1) {
                    	let gp = sim.cells.getGridpoint(ii,jj)
                        living_neighbours.push(gp) // Otherwise, add it to the array
                    }
                }
            }
            return living_neighbours
        }


		//	From referenced grid point, call findLivingNeighbours. From those gridpoints with living, return which public goods are available to this gridpoint
        sim.cells.sumPublicGoods = function (i,j) {
            let public_goods_produced = 0
            let resources = []
            resources.length = num_CG
            for (let c = 0; c < resources.length; c++) {
                resources[c] = 0
            }

            let neighbours = this.findLivingNeighbours(i,j) 
            if (neighbours.length > 0) {
                for (let z = 0; z < neighbours.length; z++) {
                    let gen_array = neighbours[z].genome
                    for (let y = 0; y < gen_array.length; y++) {
                        if (gen_array[y] == 1) {
                            resources[y] = 1
                        }
                    }
                }
            }
    
            for (let xx = 0; xx < resources.length; xx++){
                if (resources[xx] == 1) {
                    public_goods_produced += 1
                }
            }
            return public_goods_produced
        }


		//	Check if CG resources are available, roll dice for reproduction (with penalty for amount of CGs produced)
        sim.cells.checkGrowth = function(i,j) {
            if (this.grid[i][j].alive == 0) {
                let sum_resources = sim.cells.sumPublicGoods(i,j)
                if (sum_resources == num_CG) {
                    let neighbours = this.getMoore8(this, i, j, 'alive', 1)
                    let random_neighbours = neighbours
                    shuffle(random_neighbours)
                    for (let bb = 0; bb < random_neighbours.length; bb++) {
                    	let growth_penalty = cost_per_CG * base_repro_chance
                        let repro_chance = base_repro_chance - (random_neighbours[bb].production * growth_penalty)
                        if (this.rng.random() < repro_chance) {
                            
                            this.grid[i][j].alive = 1
                            this.grid[i][j].genome = new Array(num_CG)
                            for(let p = 0; p < num_CG; p++) {
                            	this.grid[i][j].genome[p] = random_neighbours[bb].genome[p]  
							}
                            temp_var = this.grid[i][j].genome
                            let prod_count = 0
                            for(let u = 0; u < temp_var.length; u++) {
                                if(temp_var[u] == 1) {
                                    prod_count += 1
                                	}                               
                        	}
                            this.grid[i][j].production = prod_count
                        }
                    }
                }
            }
        }


		//	Func used to count production for each number of CGs (i.e., number 1x producers, 2x producers, etc.) for writing out
    	sim.cells.getProdAmounts = function(property, values) {
        	let sum = Array(values.length).fill(0);
	  		for (let i = 0; i < this.nc; i++) {
    		    for (let j = 0; j < this.nr; j++) {
            		for (let val in values)
                	if (this.grid[i][j].alive == 1 ){
                		if (this.grid[i][j][property] == values[val]) sum[val]++;
                	}
	           	}
        	}
	        return sum;
    	}


		//	For the doppleganger grid; feed a reference gridpoint. Return all living cells within a box within the radius of CG effect in this sim.
        sim.dopple.findLivingNeighbours = function (i,j) {
            let living_neighbours = [] // empty array where you will store grid points in
            for(let ii= (i - CG_radius); ii<= (i + CG_radius); ii++)
            {
                for(let jj= (j - CG_radius); jj<= (j + CG_radius); jj++)
                {
                    if(sim.dopple.getGridpoint(ii,jj) == undefined) {}
                    else if(sim.dopple.getGridpoint(ii,jj).alive == 1) {
	                    let gp = sim.dopple.getGridpoint(ii,jj)
                        living_neighbours.push(gp) // Otherwise, add it to the array
                    }
                }
            }
            return living_neighbours
        }

	
		//	For the doppleganger grid; From referenced grid point, call findLivingNeighbours. From those gridpoints with living, return which public goods are available to this gridpoint
        sim.dopple.sumPublicGoods = function (i,j) {
            let public_goods_produced = 0
            let resources = []
            resources.length = num_CG
            for (let c = 0; c < resources.length; c++) {
                resources[c] = 0
            }

            let neighbours = this.findLivingNeighbours(i,j) // All neighbours in a 10 (square) radius, if you want a circular neighbourhood than you need to first calculate the distance a la pytagoras
            if (neighbours.length > 0) {
                for (let z = 0; z < neighbours.length; z++) {
                    let gen_array = neighbours[z].genome
                    for (let y = 0; y < gen_array.length; y++) {
                        if (gen_array[y] == 1) {
                            resources[y] = 1
                        }
                    }
                }
            }
    
            for (let xx = 0; xx < resources.length; xx++){
                if (resources[xx] == 1) {
                    public_goods_produced += 1
                }
            }
            return public_goods_produced
        }


		//	For the doppleganger grid; Check if CG resources are available, roll dice for reproduction (with penalty for amount of CGs produced)
        sim.dopple.checkGrowth = function(i,j) {
            if (this.grid[i][j].alive == 0) {
                let sum_resources = sim.dopple.sumPublicGoods(i,j)
                if (sum_resources == num_CG) {
                    let neighbours = this.getMoore8(this, i, j, 'alive', 1)
                    let random_neighbours = neighbours
                    shuffle(random_neighbours)
                    for (let bb = 0; bb < random_neighbours.length; bb++) {
                    	let growth_penalty = cost_per_CG * base_repro_chance
                        let repro_chance = base_repro_chance - (random_neighbours[bb].production * growth_penalty)
                        if (this.rng.random() < repro_chance) {
                            this.grid[i][j].alive = 1
                            this.grid[i][j].genome = new Array(num_CG)
                            for(let p = 0; p< num_CG; p++) {
                            	this.grid[i][j].genome[p] = random_neighbours[bb].genome[p]  
							}
                            temp_var = this.grid[i][j].genome
                            let prod_count = 0
                            for(let u = 0; u < temp_var.length; u++) {
                                if(temp_var[u] == 1) {
                                    prod_count += 1
                                }
                            }
                            this.grid[i][j].production = prod_count
                        }
                    }
                }
            }
        }


		//	For the doppleganger grid; Func used to count production for each number of CGs (i.e., number 1x producers, 2x producers, etc.) for writing out
    	sim.dopple.getProdAmounts = function(property, values) {
        	let sum = Array(values.length).fill(0);
	  		for (let i = 0; i < this.nc; i++) {
    		    for (let j = 0; j < this.nr; j++) {
            		for (let val in values)
                	if (this.grid[i][j].alive == 1 ){
                		if (this.grid[i][j][property] == values[val]) sum[val]++;
                	}
	           	}
        	}
	        return sum;
    	}

    	
    	//	function to initiate invasion. Swaps everything in a region of the cells grid with the same region of the dopple grid. Creates a double invasion experiment by continuing both; currently set up as a 50/50 swap.
    	sim.cells.invasion_byregion = function(i,j) {
    		for (let i = 0; i < sim.cells.nc; i++) for (let j = 0; j < sim.cells.nr; j++) {			
	            if ((i < (sim.cells.nc / 2)) && (j < (sim.cells.nr))) {
                    temp_alive = sim.cells.grid[i][j].alive
                    temp_genome = sim.cells.grid[i][j].genome
                    temp_prod = sim.cells.grid[i][j].production
					
                    sim.cells.grid[i][j].alive = sim.dopple.grid[i][j].alive 
                    sim.cells.grid[i][j].genome = sim.dopple.grid[i][j].genome
                    sim.cells.grid[i][j].production = sim.dopple.grid[i][j].production
					
                    sim.dopple.grid[i][j].alive = temp_alive 
                    sim.dopple.grid[i][j].genome = temp_genome
                    sim.dopple.grid[i][j].production = temp_prod              	                    
        	    }
	        }
    	}        			
			
    	

// Write out functions
	// Write out Headers for each cheater/producer permutation


// Doppleganger write out functions

		sim.dopple.writeFirstCheaterPermutations = function() {
			let filn = dir+"/CheaterCounts_"+"doppleganger_"+filn_suffix+".txt"	
			fs.appendFileSync(filn, `TimeSteps\t`)	
			for (let [key, value] of masterPerm) {
				fs.appendFileSync(filn, `${key}\t`)
			}
			fs.appendFileSync(filn, `\n`)
		}
	// Write out counts for each cheater/producer permutation
		sim.dopple.writeCheaterPermutations = function() {
			let filn = dir+"/CheaterCounts_"+"doppleganger_"+filn_suffix+".txt"		
			fs.appendFileSync(filn, `${this.time}\t`)	
			for (let [key, value] of permCounts) {
				fs.appendFileSync(filn, `${value}\t`)
			}
			fs.appendFileSync(filn, `\n`)	
		}

	// Write grid
		sim.dopple.writeGrid = function() {
			filn = dir+"/GridStates_"+"doppleganger_"+filn_suffix+"/"+this.time+".txt"	
			if(fs.existsSync(filn)){
				fs.unlinkSync(filn)
			}
			fs.appendFileSync(filn, "x\ty\talive\tproduction_count\tgenome_state\n")
			for (let i = 0; i < this.nc; i++) {
    		    for (let j = 0; j < this.nr; j++) {    		    	
    		    	let point_state = i+'\t'+j+'\t'+this.grid[i][j].alive+'\t'+this.grid[i][j].production+'\t'+this.grid[i][j].genome+'\n'
 					fs.appendFileSync(filn,point_state)	
 				}
 			}	
		}


// Normal write out functions (cells)
		sim.cells.writeFirstCheaterPermutations = function() {
			let filn = dir+"/CheaterCounts_"+filn_suffix+".txt"	
			fs.appendFileSync(filn, `TimeSteps\t`)	
			for (let [key, value] of masterPerm) {
				fs.appendFileSync(filn, `${key}\t`)
			}
			fs.appendFileSync(filn, `\n`)
		}
	// Write out counts for each cheater/producer permutation
		sim.cells.writeCheaterPermutations = function() {
			let filn = dir+"/CheaterCounts_"+filn_suffix+".txt"	
			fs.appendFileSync(filn, `${this.time}\t`)	
			for (let [key, value] of permCounts) {
				fs.appendFileSync(filn, `${value}\t`)
			}
			fs.appendFileSync(filn, `\n`)	
		}

	// Write out Headers for alive cells count
		sim.cells.writeFirstAlive = function() {
			let filn = dir+"/count_alive_"+filn_suffix+".txt"
			fs.appendFileSync(filn, `TimeSteps\tAlive_count\n`)				
		
		}
	// Write out count of alive cells
		sim.cells.writeAlive = function() {
			let filn = dir+"/count_alive_"+filn_suffix+".txt"
						
			let out_2 = count_alive+'\n'
			let out_1 = this.time+'\t'
			fs.appendFileSync(filn, out_1)
			fs.appendFileSync(filn, out_2)		
		}

	// Write out Headers for counts of different CG producers
		sim.cells.writeFirstProd = function() {
			let filn = dir+"/production_"+filn_suffix+".txt"		
			let prod_arr = []
			for(let arr = 0; arr <= num_CG; arr++) {
				prod_arr[arr] = arr
			}
			fs.appendFileSync(filn, 'TimeSteps\t')
			for(let p = 0; p <= num_CG; p++) {
				fs.appendFileSync(filn, "Total_CGs_"+p+'\t')
			}
			fs.appendFileSync(filn, '\n')
		}
	// Write out counts of different CG producers
		sim.cells.writeProd = function() {
			let filn = dir+"/production_"+filn_suffix+".txt"		
			let prod_arr = []
			for(let arr = 0; arr <= num_CG; arr++) {
				prod_arr[arr] = arr
			}
			let production_counts = sim.cells.getProdAmounts('production', prod_arr)
			fs.appendFileSync(filn, this.time+'\t')
			for(let p = 0; p <= num_CG; p++) {
				fs.appendFileSync(filn, production_counts[p]+'\t')
			}
			fs.appendFileSync(filn, '\n')
		}
	// Write Headers for number of deaths
		sim.cells.writeFirstKillCount = function() {
			let filn = dir+"/kill_count_"+filn_suffix+".txt"						
			let out_2 = 'Kills\n'
			let out_1 = 'TimeSteps\t'
			fs.appendFileSync(filn, out_1)
			fs.appendFileSync(filn, out_2)		
		}	
	// Write number of deaths	
		sim.cells.writeKillCount = function() {
			let filn = dir+"/kill_count_"+filn_suffix+".txt"					
			let out_2 = kills+'\n'
			let out_1 = this.time+'\t'
			fs.appendFileSync(filn, out_1)
			fs.appendFileSync(filn, out_2)		
		}				
	// Write grid
		sim.cells.writeGrid = function() {
			filn = dir+"/GridStates_"+filn_suffix+"/"+this.time+".txt"	
			if(fs.existsSync(filn)){
				fs.unlinkSync(filn)
			}
			fs.appendFileSync(filn, "x\ty\talive\tproduction_count\tgenome_state\n")
			for (let i = 0; i < this.nc; i++) {
    		    for (let j = 0; j < this.nr; j++) {    		    	
    		    	let point_state = i+'\t'+j+'\t'+this.grid[i][j].alive+'\t'+this.grid[i][j].production+'\t'+this.grid[i][j].genome+'\n'
 					fs.appendFileSync(filn,point_state)	
 				}
 			}	
		}

        sim.cells.initialise()
        sim.cells.writeGrid()
		First_makePermCounts()		// Create the first PermCount Hash that will guide printing of future hash updates   
		sim.cells.writeFirstProd()
		sim.cells.writeFirstCheaterPermutations()
        sim.dopple.initialise()
        sim.dopple.writeGrid()
        sim.dopple.writeFirstCheaterPermutations()
        sim.dopple.writeCheaterPermutations()
        sim.cells.writeCheaterPermutations()


		//	Main grid's next state function
        sim.cells.nextState = function (i, j) {

            if (this.grid[i][j].alive == 1) {

                // Mutate at global var's rate 
                for (let t = 0; t < num_CG; t++) {
                    let rng = this.rng.random()
                    if ((rng < global_mut) && (this.grid[i][j].genome[t] == 1)) {
                        this.grid[i][j].genome[t] = 0
                    	}                    
				}

				// HGT from neighbours using Chromosome_HGT_Rate
                for (let n_CG = 0; n_CG < num_CG; n_CG++) {
                    let rng = this.rng.random()
					if (((rng < (Chromosome_HGT_rate * Zero_HGT_ratio)) && (this.grid[i][j].genome[n_CG] == 0)) || ((rng < Chromosome_HGT_rate) && (this.grid[i][j].genome[n_CG] == 1))) {	
						let neighbours = this.getMoore8(this, i, j, 'alive', 1)
						if (neighbours.length >= 1) {
							shuffle(neighbours)
							neighbours[0].genome[n_CG] = this.grid[i][j].genome[n_CG]
						}
					}
				}

                temp_var = this.grid[i][j].genome
                this.grid[i][j].production = temp_var.reduce((a, b) => a + b, 0)

                // Kill at global var's rate
                if (this.rng.random() < global_death) {
                    kills += 1
                    this.grid[i][j].alive = 0
                    this.grid[i][j].genome = undefined
                    this.grid[i][j].production = 0
                }
                else count_alive ++
            }

            // Check for growth
            sim.cells.checkGrowth(i,j)
		}


		// Doppleganger grid's NextState
        sim.dopple.nextState = function (i, j) {

            if (this.grid[i][j].alive == 1) {
                // Mutate at global var's rate 
                for (let t = 0; t < num_CG; t++) {
                    let rng = this.rng.random()
                    if ((rng < global_mut) && (this.grid[i][j].genome[t] == 1)) {
                        this.grid[i][j].genome[t] = 0
                    	}                
				}
 
 				// HGT from neighbours using Chromosome_HGT_Rate
                for (let n_CG = 0; n_CG < num_CG; n_CG++) {
                    let rng = this.rng.random()
					if (((rng < (Chromosome_HGT_rate * Zero_HGT_ratio)) && (this.grid[i][j].genome[n_CG] == 0)) || ((rng < Chromosome_HGT_rate) && (this.grid[i][j].genome[n_CG] == 1))) {	
						let neighbours = this.getMoore8(this, i, j, 'alive', 1)
						if (neighbours.length >= 1) {
							shuffle(neighbours)
							neighbours[0].genome[n_CG] = this.grid[i][j].genome[n_CG]
						}
					}
				}

                temp_var = this.grid[i][j].genome
                this.grid[i][j].production = temp_var.reduce((a, b) => a + b, 0)

                // Kill at global var's rate
                if (this.rng.random() < global_death) {
                    kills += 1
                    this.grid[i][j].alive = 0
                    this.grid[i][j].genome = undefined
                    this.grid[i][j].production = 0
                }
                
            }
            // Check for growth
            sim.dopple.checkGrowth(i,j)
        }



        /*
            3. MAIN SIMULATION LOOP. 
        */

        sim.cells.update = function () {
            this.synchronous()         // Applied as many times as it can in 1/60th of a second
            if(sim.time % 50 == 0) sim.log(`At T=${sim.time} there are ${count_alive} individuals alive`,"output")
            if (diffuse > 0){
            	this.MargolusDiffusion()
            }
 
            if (this.time % write_out_interval == 0)       // Time zero print out
            {
				sim.cells.writeProd()
				//Update CheaterPermutations, and then write
				sim.cells.reportGenomes()
				sim.cells.writeCheaterPermutations()
            }

            if (this.time == 15000)
            {
//            	sim.cells.invasion_byregion()
            }

            if (this.time % (write_out_interval * 5) == 0)		// write out grid state less frequently than reporting genotypes
            {
            	sim.cells.writeGrid()
            } 
                      
			count_alive = 0
        }
        
        
        sim.dopple.update = function () {
            this.synchronous()         // Applied as many times as it can in 1/60th of a second
            if (diffuse > 0){
            	this.MargolusDiffusion()
            }

            if (this.time % write_out_interval == 0)       // Time zero print out
            {
				//Update CheaterPermutations, and then write
				sim.dopple.reportGenomes()
				sim.dopple.writeCheaterPermutations()		
            }
            
            if (this.time == 5000)			// trigger invasion by region
            {
//            	sim.dopple.invasion_byregion()
            }

            if (this.time % (write_out_interval * 5) == 0)		// write out grid state less frequently than reporting genotypes
            {
            	sim.dopple.writeGrid()
            } 
                    
			count_alive = 0
        }      

        sim.start()
    }
        
        
if (typeof window == "undefined") cacatoo()