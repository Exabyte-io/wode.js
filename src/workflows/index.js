import { allApplications } from "@exabyte-io/ade.js";

import { Workflow } from "./workflow";
import { createWorkflow } from "./create";
import { workflowData as allWorkflowData } from "./workflows";


/*
    Workflow construction follows these rules:
        1. Workflow is constructed as a collection of subworkflows defined in JSON
        2. A "units" key should contain at least one object referencing the workflow itself
        3. Additional workflows are added in order specified in the same "units" array
        4. map units are added along with their workflows according to data in "units"
        5. top-level subworkflows are added directly in the order also specified by "units"
 */
function createWorkflows({
    appName = null,
    workflowCls = Workflow,
    ...swArgs
}) {
    const apps = (appName !== null) ? [appName] : allApplications;
    const wfs = [];
    const { workflows } = allWorkflowData;
    for (const appName of apps) {
        const { [appName]: dataByApp } = workflows;
        for (const workflowData of Object.values(dataByApp)) {
            wfs.push(createWorkflow({
                appName, workflowData, workflowCls, ...swArgs,
            }));
        }
    }
    return wfs;
}

export { Workflow, createWorkflows };
