// Tool class definition
class Tool {
    String name
    List<String> compatibleTasks
	List<String> compatibleDatasets

    Tool(String name, List<String> compatibleTasks, List<String> compatibleDatasets) {
        this.name 				= name
        this.compatibleTasks 	= compatibleTasks
		this.compatibleDatasets = compatibleDatasets
    }

    // Check if the tool is compatible with a given perturbation task
    boolean isTaskCompatible(String task) {
        return this.compatibleTasks.contains(task)
    }

	// Check if the tool is compatible with a given dataset type
    boolean isDatasetCompatible(String dataset_type) {
        return this.compatibleDatasets.contains(dataset_type)
    }

    String toString() {
        "$name"
    }
}