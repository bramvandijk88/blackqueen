      
let num_CG = 6			// Number of CGs in system
var global_mut = 0.00005		//	Chance for each gene in each cell to mutate into non-production per timestep
let global_death = 0.1		// Chance of death for each cell per timestep
let base_repro_chance = 0.6		// Default chance a non-producer may be able to reproduce
let starting_omni_proportion = 0.6
var cost_per_CG = 0.0001		// Proportionate penalty for each CG 

var CG_radius = 1		// radius of Moore neighbourhood for benefits from common goods
let temp_var
let count_alive = 0			// How many cells are alive
let births = 0						// Count num births
let kills = 0 
let maxsteps = 1000000                             // How many time steps the model continues to run            
let num_col = 160                                   // Number of columns (width of your grid)
let num_row = 160								// Number of columns (width of your grid)
let diffuse = 1				// If >0, then margolis diffusion will be enabled

var Chromosome_HGT_rate = 0		// Chance of a CG gene in a living cell to copy itself over its ortholog in a neighbouring cell's genome
let Zero_HGT_ratio = 1		// Bias in HGT of non-producing alleles. 1 = zeros transfer at same rate as 1 (producing gene); 0 = zeros never transfer (treats them as deleted over inactivated)

var list = [0,1]
let permCounts = undefined
let masterPerm = undefined
let run_num = 1

var display_interval = 0
var mix = 0


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
        title: "BQ dynamics",                        // The name of your cacatoo-simulation
        description: "",         // description, if you needed
        maxtime: maxsteps,                             // How many time steps the model continues to run            
        ncol: num_col,                                   // Number of columns (width of your grid)
        nrow: num_row,		                              // Number of rows (height of your grid)
        graph_interval: 20,
        graph_update: 100,
        seed: 2,
        skip:display_interval,
        wrap: [true, true],                         // Wrapped boundary conditions? [COLS, ROWS]   
        scale: 2,				                      // Scale of the grid (nxn pixels per grid point)
        bgcolour: '#000000'
    }

    /*
        1. SETUP. (continued) 
    */
    sim = new Simulation(config)                          // Initialise the Cacatoo simulation
    sim.makeGridmodel("cells")                              // Build a new Gridmodel within the simulation called "cells"
    sim.makeGridmodel("frequencies")
    sim.frequencies.nc = num_CG+1
    sim.frequencies.nr = num_CG
    sim.frequencies.scale = Math.round(config.scale*config.ncol/num_CG)

    sim.cells.colourGradient('production', num_CG+1, [0, 30, 0], [0, 200, 0])                         // Will contain the ODEs, and show the abundance of PREDATORS 
    sim.cells.colourGradient('available', num_CG+1, [0, 30, 0], [0, 200, 0])                         // Will contain the ODEs, and show the abundance of PREDATORS 
    sim.frequencies.colourGradient('value', num_CG+1, [0, 30, 0], [0, 200, 0])                         // Will contain the ODEs, and show the abundance of PREDATORS 

    //sim.createDisplay_discrete({model:"cells", property:"alive", label:"Alive"})
    
    sim.createDisplay_continuous({model:"cells", property:"production", label:"#CGs <b>produced</b> (cells)", 
                                    minval:0, maxval:num_CG, decimals: 0, nticks:2, num_colours:num_CG+1
                                    })
    sim.createDisplay_continuous({model:"cells", property:"available", label:"#CGs <b>available</b> (environment)",
                                    minval:0, maxval:num_CG, decimals: 0, nticks:2, num_colours:num_CG+1,drawdots:false})
    
    sim.createDisplay_continuous({model:"frequencies", property:"value", label:`Ecosystem table (top ${num_CG+1} shown)`,
                        minval:0, maxval:num_CG, decimals: 0, nticks:2, num_colours:num_CG+1,drawdots:false})
    
    /*
        2. DEFINING THE RULES. Below, the user defines the nextState function. This function will be applied for each grid point when we will update the grid later. 
    */

    //	Initialise a grid - write an initial accounting through reportGenomes() function
    sim.cells.initialise = function () {
        sim.initialGrid(sim.cells, "production", 0, 1.0)
        sim.initialGrid(sim.cells, "available", 1, 1.0)
        sim.initialGrid(sim.frequencies, "value", 0, 1.0)
        for (let i = 0; i < sim.cells.nc; i++) for (let j = 0; j < sim.cells.nr; j++) {
            if (sim.cells.rng.random() < starting_omni_proportion) {
                sim.cells.grid[i][j].alive = 1
                sim.cells.grid[i][j].genome = new Array(num_CG).fill(1)             
                temp_var = sim.cells.grid[i][j].genome
                sim.cells.grid[i][j].production = temp_var.reduce((a, b) => a + b, 0)
            }
           else sim.cells.grid[i][j].alive=0
        }
        sim.cells.reportGenomes()
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
        permCounts = new Map([...permCounts.entries()].sort((a, b) => b[1] - a[1]));
    }   

    
          
    
    
    //	feed a reference gridpoint. Return all living cells within a box within the radius of CG effect in this sim.
    sim.cells.findLivingNeighbours = function (i,j) {
        let living_neighbours = [] // empty array where you will store grid points in
        for(let ii= (i - CG_radius); ii<= (i + CG_radius); ii++)
        {
            for(let jj= (j - CG_radius); jj<= (j + CG_radius); jj++)
            {
                if(sim.cells.getGridpoint(ii,jj).alive == 1) {
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
        this.grid[i][j].available = public_goods_produced
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
     

    sim.cells.initialise()

    //	Main grid's next state function
    sim.cells.nextState = function (i, j) {
        
        if (this.grid[i][j].alive == 1) {
            count_alive++
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
                this.grid[i][j].production = undefined
                this.grid[i][j].available = 0
            }            
        }

        // Check for growth
        sim.cells.checkGrowth(i,j)
    }


    

    /*
        3. MAIN SIMULATION LOOP. 
    */
    sim.frequencies.update = function(){
         
    }
    sim.frequencies.heatmap = function(){
        let n=0
        for (var m of permCounts){
            var values = m[0].split(',').map(function(item) {
                return parseInt(item, 10);
            });
            let sum = values.reduce((a, b) => a + b, 0)
            
            for(let i=0; i<num_CG+1; i+=1)
                if(values[i]>0) this.grid[n][i].value = sum
                else this.grid[n][i].value = 0
            n++
            if(n == num_CG+1) return;
        }  
    }

    sim.cells.update = function () {
        this.synchronous()         // Applied as many times as it can in 1/60th of a second
        if(mix) this.perfectMix()
        if (diffuse > 0){
            this.MargolusDiffusion()
        }

        if (this.time == 15000)
        {
//            	sim.cells.invasion_byregion()
        }

        // Output stuff
        let sumproduction_array = []
        let production_array = []
        let production_cols = []
        let production_labs = []
        for(let i=0;i<num_CG+1;i++) {
            production_labs[i] = `${i} producer`
            production_array[i] = 0
            production_cols[i] = [0,Math.round(200*i/num_CG),0]
        }
        for (let x = 0; x < this.nc; x++)          // x are columns
            for (let y = 0; y < this.nr; y++)      // y are rows
            {
                if(this.grid[x][y].alive == 1) {
                    sumproduction_array.push(this.grid[x][y].production)
                    production_array[this.grid[x][y].production] += 1
                }
            }
        sumproduction_array = shuffle(sumproduction_array)
        sumproduction_array = sumproduction_array.slice(-100)
       // this.plotPoints(sumproduction_array,"CG production in a sample of 100 individuals",{width:800, height:300, labelsDivWidth: 0})   
        this.plotArray(production_labs,production_array,production_cols,"Frequencies of n-producers",{width:800, height:300, labelsDivWidth: 200})   
        
        if(sim.time%100==0){
            sim.frequencies.heatmap()
            sim.cells.reportGenomes()
        }
             
        count_alive = 0
    }
    
    
    
    //sim.addButton("step", function () { sim.step(); sim.display() })  
    
    sim.addSlider("global_mut",0,0.01,0.00001,"Mutation rate")
    sim.addSlider("CG_radius",1,5,1,"Interaction range")
    sim.addSlider("cost_per_CG",0,0.1,0.0001,"CG production costs")
    sim.addSlider("Chromosome_HGT_rate",0,0.002,0.00001,"HGT rate")
    sim.addHTML("form_holder","<br>")
    sim.addButton("Pause/continue", function () { sim.toggle_play() })
    sim.addButton("Download data", function() { 
        let data = sim.cells.graphs["Frequencies of n-producers"].data
        let str = 'Time'
        for(let i=0;i<=num_CG;i++) str+= `,${i}-producer`
        str += '\n'
        for(let i=0;i<data.length; i++)
            if(data[i][0]%100==0)  str+= data[i] + '\n'
        sim.write(str,"Timeseries_BQ.txt")
    }) 
    //sim.addCustomSlider("Display interval",function (new_value) { sim.skip = new_value },0,100,5,sim.skip)
    
    sim.start()
}
        
        
if (typeof window == "undefined") cacatoo()
