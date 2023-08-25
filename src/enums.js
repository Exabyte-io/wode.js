/**
 * THIS ENUMS ARE SHARED WITH TESTS.
 * DO NOT IMPORT ANYTHINGS IN THIS MODULE.
 */

export const IO_ID_COLUMN = "exabyteId";

export const UNIT_TYPES = {
    // not currently used
    convergence: "convergence",
    exit: "exit",
    // actively used
    execution: "execution",
    map: "map",
    reduce: "reduce",
    assignment: "assignment",
    condition: "condition",
    subworkflow: "subworkflow",
    processing: "processing",
    io: "io",
    assertion: "assertion",
};

export const UNIT_STATUSES = {
    idle: "idle",
    active: "active",
    finished: "finished",
    error: "error",
    warning: "warning",
};

export const UNIT_TAGS = {
    hasConvergenceParam: "hasConvergenceParam",
    hasConvergenceResult: "hasConvergenceResult",
};

export const WORKFLOW_STATUSES = {
    "up-to-date": "up-to-date",
    outdated: "outdated",
};

export const TAB_NAVIGATION_CONFIG = {
    overview: {
        itemName: "Overview",
        className: "",
        href: "sw-overview",
    },
    importantSettings: {
        itemName: "Important settings",
        className: "",
        href: "sw-important-settings",
    },
    detailedView: {
        itemName: "Detailed view",
        className: "",
        href: "sw-detailed-view",
    },
    compute: {
        itemName: "Compute",
        className: "",
        href: "sw-compute",
    },
};

export const UNIT_NAME_INVALID_CHARS = "/";
