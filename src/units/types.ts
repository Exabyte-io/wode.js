export enum UNIT_TYPES {
    execution = "execution",
    map = "map",
    reduce = "reduce",
    assignment = "assignment",
    condition = "condition",
    subworkflow = "subworkflow",
    processing = "processing",
    io = "io",
    assertion = "assertion",
};

export enum UNIT_STATUSES {
    idle = "idle",
    active = "active",
    finished = "finished",
    error = "error",
    warning = "warning",
};

export type UnitInput = object[];

export interface NamedEntity {
    name: string;
}

export interface TaggableEntity {
    tags: string[];
}

export interface StatusTrack {
    repetition: number;
    trackedAt: string; // ISO date string
}

export interface StatusTrackingEntity {
    status: string;
    statusTrack: StatusTrack[];
}

export interface UnitConfig extends NamedEntity, TaggableEntity, StatusTrackingEntity {
    flowchartId: string;
    type: UNIT_TYPES;
    head: boolean; // true if the unit is the head of the flowchart
    next: string;
    isDraft: boolean; // TODO: figure out where this is used
}

export interface AssertionUnitConfig extends UnitConfig {
    type: UNIT_TYPES.assertion;
    statement: string;
    errorMessage: string;
}

export interface AssignmentUnitConfig extends UnitConfig {
    type: UNIT_TYPES.assignment;
    input: UnitInput;
    operand: string;
    value: string;
}

export interface ConditionUnitConfig extends UnitConfig {
    type: UNIT_TYPES.condition;
    statement: string;
    maxOccurrences: number;
    then: string;
    else: string;
}

// TODO: define Executable and Flavor types
export interface ApplicationType extends NamedEntity {
    defaultExecutable: ExecutableType;
    getExecutableByConfig: (config: ExecutableConfig) => ExecutableType;
    toJSON: () => this
};

export interface ExecutableConfig {any};
export interface ExecutableType extends NamedEntity {
    defaultFlavor: FlavorType;
    results: string[]; // array of result names
    monitors: string[]; // array of monitor names
    preProcessors: string[]; // array of preprocessor names
    postProcessors: string[] // array of postprocessor names
    getFlavorByConfig: (config: FlavorConfig) => FlavorType;
    toJSON: () => this
};

export interface TemplateType {
    getContextProvidersAsClassInstances: (contextProviders: any[]) => any[];
};
export interface FlavorConfig {any};
export interface FlavorType extends NamedEntity {
    inputAsTemplates: TemplateType[];
    results: string[]; // array of result names
    monitors: string[]; // array of monitor names
    preProcessors: string[]; // array of preprocessor names
    postProcessors: string[] // array of postprocessor names
    getInputAsRenderedTemplates: (context: any) => TemplateType[];
    toJSON: () => this
};

export interface ExecutionUnitConfig extends UnitConfig {
    type: UNIT_TYPES.execution;
    application: ApplicationType;
    executable: ExecutableType;
    flavor: FlavorType;
    input: UnitInput;
}

export type IOUnitSubTypes = 'input' | 'output' | 'dataFrame';
export interface IOUnitConfig extends UnitConfig {
    type: UNIT_TYPES.io;
    subtype: IOUnitSubTypes;
    source: string;
    input: UnitInput;
    enableRender: boolean; // whether template should be rendered with jinja templates at runtime
}

export interface MapUnitConfig extends UnitConfig {
    type: UNIT_TYPES.map;
    input: UnitInput;
    workflowId: string;
}

export interface ProcessingUnitConfig extends UnitConfig {
    type: UNIT_TYPES.processing;
    operation: string;
    operationType: string;
    input: UnitInput;
}

export interface ReduceUnitConfig extends UnitConfig {
    type: UNIT_TYPES.reduce;
    mapFlowchartId: string;
    input: UnitInput;
}

export interface SubworkflowUnitConfig extends UnitConfig {
    type: UNIT_TYPES.subworkflow;
}
