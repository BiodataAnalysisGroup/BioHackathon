// Run class definition
class Run {
    Tool tool
    Dataset dataset
    Task task

    Run(Tool tool, Dataset dataset, Task task){
        this.tool       = tool
        this.dataset    = dataset
        this.task       = task
    }

    String toString() {
        "Run for tool $tool on dataset $dataset for task $task"
    }

    String getToolName(){
        "$tool"
    }
}