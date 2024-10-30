#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2


/*
================================================================
	MODULES
================================================================
*/
// Include modules
include {PREPROCESS_DATASET} from "./modules/preprocess_dataset"
include {CELL_ORACLE} from "./modules/tools/cell_oracle"
include {CELL_OT} from "./modules/tools/cell_ot"
include {ACCURACY} from "./modules/metrics/accuracy"
include {COSINE} from "./modules/metrics/cosine"
include {SUMMARIZE_RESULTS} from "./modules/summarize_results"



/*
================================================================
	SUBWORKFLOWS
================================================================
*/
//TODO: include subworkflows




/*
================================================================
DEFINE GLOBAL VARIABLES
================================================================
*/
// Define available metrics
def all_metrics = params.metrics.collect { params -> 
    new Metric(*params)
}

// Define available datasets
def all_datasets = params.datasets.collect { params ->
    new Dataset(*params) 
}

// Define available tools
def all_tools = params.tools.collect { params ->
    new Tool(*params)
}


// Define available Task with their compatible metrics
def all_task = params.tasks.collect { params ->
    new Task(*params)
}



/*
================================================================
MAIN WORKFLOW
================================================================
*/
workflow {

	//**********
	//*  INIT  *
	//**********
    // Get selected tools or default to all
    def selected_tools = params.selected_tools ?
        all_tools.findAll { tool -> tool.name in params.selected_tools } :
        all_tools


    // Get selected tasks or default to compatible tasks
    def compatibleTasks = all_tools.collect { it.compatibleTasks }.flatten().unique()     
    def selected_tasks  = params.selected_tasks ?
        all_task.findAll { it.name in params.selected_tasks } :
        all_task.findAll { it.name in compatibleTasks } 


    // Get selected datasets or default to compatible datasets
    def compatibleDatasets = all_tools.collect { it.compatibleDatasets }.flatten().unique()
    def selected_datasets = params.selected_datasets ?
        all_datasets.findAll { it.name in params.selected_datasets } :
        all_datasets.findAll { it.type in compatibleDatasets }
    

    // Get selected metrics or default to all
    def compatibleMetrics = all_task.collect { it.compatibleMetrics }.flatten().unique()
    def selected_metrics = params.selected_metrics ?
        all_metrics.findAll { it.name in params.selected_metrics } :
        all_metrics.findAll { it.name in compatibleMetrics }
    

	// init channels
	selected_tools_ch 		= Channel.from(selected_tools)
	selected_tasks_ch 		= Channel.from(selected_tasks)
	selected_datasets_ch	= Channel.from(selected_datasets)
    selected_metrics_ch     = Channel.from(selected_metrics)



	//**********
	//*  MAIN  *
	//**********

    // Preprocess datasets (download, unzip, reformat, ...)
	processed_datasets_ch = selected_datasets_ch
        .ifEmpty { error "No datasets selected for execution." } 
        .dump(tag: 'PREPROCESSING DATABASET: ')        
        .map { dataset -> [dataset, dataset.path] } 
        | PREPROCESS_DATASET


	// Prepare running experiments
    selected_tools_ch                
        .combine(selected_tasks_ch)        
        .combine(processed_datasets_ch)         
        .filter { tool, task, dataset, dataset_path ->    
            tool.isDatasetCompatible(dataset.type) 
            tool.isTaskCompatible(task.name) }
        .tap { tool_launcher_ch } 
        .combine(selected_metrics_ch)
        .filter { tool, task, dataset, dataset_path, metric ->              
            task.isMetricCompatible(metric.name)}                 
        .tap { metric_launcher_ch }

        
    // create channel for running tools
    runs_ch = tool_launcher_ch
        .map { tool, task, dataset, dataset_path ->            
            [tool.name, task.name, dataset.name, dataset_path] }
        .dump (tag: 'RUN TOOL: ')


    // Run tools for perturbation tasks on specific datasets
    tool_outputs = CELL_ORACLE(runs_ch)
        .concat(CELL_OT(runs_ch))
        //.concat(TOOL_X(runs_ch))        
        .dump(tag: 'TOOL OUTPUT: ')


    // create channel for running metrics
    metrics_ch = metric_launcher_ch		
        .map {tool, task, dataset, dataset_path, metric ->               
            [tool.name, task.name, dataset.name, metric.name]}
		.combine( tool_outputs, by: [0,1,2] )  // tool_name, task_name and database_name
        .dump (tag: 'RUN METRICS: ')


    // Run metrics for tool outputs on specific task/dataset
    metric_outputs = ACCURACY(metrics_ch)
		.concat(COSINE(metrics_ch))
		//.concat(METRIC_X(metrics_ch))
		.dump(tag: 'METRIC OUTPUT: ')

	
	// Summarize results
	metric_outputs
		.groupTuple(by: 1)	// group by perturbation task
		|SUMMARIZE_RESULTS



}	