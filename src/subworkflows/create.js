import ApplicationRegistry from "@exabyte-io/ade.js/dist/js/ApplicationRegistry";
import {
    default_methods as MethodConfigs,
    default_models as ModelConfigs,
    MethodFactory,
    ModelFactory,
} from "@exabyte-io/mode.js";
import _ from "lodash";

import { UnitFactory } from "../units";
import { builders } from "../units/builders";
import { applyConfig } from "../utils";
import { workflowData as allWorkflowData } from "../workflows/workflows";
import { dynamicSubworkflowsByApp, getSurfaceEnergySubworkflowUnits } from "./dynamic";
import { Subworkflow } from "./subworkflow";

// NOTE: DFTModel => DFTModelConfig, configs should have the same name as the model/method class + "Config" at the end
function _getConfigFromModelOrMethodName(name, kind) {
    const configs = kind === "Model" ? ModelConfigs : MethodConfigs;
    if (!configs[`${name}Config`]) {
        // eslint-disable-next-line no-param-reassign
        name = `Unknown${kind}`;
    }
    return configs[`${name}Config`];
}

/**
 * @summary Create model from subworkflow data
 * @param config {Object} model config
 * @param modelFactoryCls {any} model factory to use
 * @returns {DFTModel|Model}
 */
function createModel({ config, modelFactoryCls }) {
    const { name, config: modelConfig = {} } = config;
    const defaultConfig = _getConfigFromModelOrMethodName(name, "Model");
    return modelFactoryCls.create({ ...defaultConfig, ...modelConfig });
}

/**
 * @summary Create method from subworkflow data
 * @param config {Object} method configuration
 * @param methodFactoryCls {any}
 * @returns {{method, setSearchText}}
 */
function createMethod({ config, methodFactoryCls }) {
    const { name, setSearchText = null, config: methodConfig = {} } = config;
    const defaultConfig = _getConfigFromModelOrMethodName(name, "Method");
    const method = methodFactoryCls.create({ ...defaultConfig, ...methodConfig });
    return { method, setSearchText };
}

/**
 * @summary Create top-level objects used in subworkflow initialization
 * @param subworkflowData {Object} subworkflow data
 * @param AppRegistry
 * @param modelFactoryCls {any} model factory class
 * @param methodFactoryCls {any} method factory class
 * @returns {{application: *, method: *, model: (DFTModel|Model), setSearchText: String|null}}
 */
function createTopLevel({ subworkflowData, modelFactoryCls, methodFactoryCls, AppRegistry }) {
    const { application: appConfig, model: modelConfig, method: methodConfig } = subworkflowData;
    const application = AppRegistry.createApplication(appConfig);
    const model = createModel({ config: modelConfig, modelFactoryCls });
    const { method, setSearchText } = createMethod({ config: methodConfig, methodFactoryCls });
    return {
        application,
        model,
        method,
        setSearchText,
    };
}

/**
 * @summary Create workflow unit from JSON configuration
 *      Supports applying functions to the builder prior to building via "functions"
 *      Supports applying attributes to the builder after building via "attributes"
 * @param config {Object} unit config
 * @param application {*} application
 * @param unitBuilders {Object} workflow unit builders
 * @param unitFactoryCls {*} workflow unit class factory
 * @returns {*|{head: boolean, preProcessors: [], postProcessors: [], name: *, flowchartId: *, type: *, results: [], monitors: []}}
 */
export function createUnit({ config, application, unitBuilders, unitFactoryCls }) {
    const { type, config: unitConfig } = config;
    if (type === "executionBuilder") {
        const { name, execName, flavorName, flowchartId } = unitConfig;
        const builder = new unitBuilders.ExecutionUnitConfigBuilder(
            name,
            application,
            execName,
            flavorName,
            flowchartId,
        );

        // config should contain "functions" and "attributes"
        const cfg = applyConfig({ obj: builder, config, callBuild: true });
        return unitFactoryCls.create(cfg);
    }

    return unitFactoryCls.create({ type, ...unitConfig });
}

/**
 * @summary Dynamically create subworkflow units
 * @param dynamicSubworkflow {String} name of unit creation function
 * @param units {Array} configured units to provide to dynamic unit creation
 * @param unitBuilders {Object} unit configuration builders
 * @param unitFactoryCls {*} unit factory class
 * @param application {*} application (optional)
 * @returns {*}
 */
function createDynamicUnits({
    dynamicSubworkflow,
    units,
    unitBuilders,
    unitFactoryCls,
    application = null,
}) {
    const { name, subfolder } = dynamicSubworkflow;
    const func = subfolder && _.get(dynamicSubworkflowsByApp, `${subfolder}.${name}`, () => {});
    switch (name) {
        case "surfaceEnergy":
            // eslint-disable-next-line no-case-declarations
            const [scfUnit] = units;
            return getSurfaceEnergySubworkflowUnits({ scfUnit, unitBuilders });
        case "getQpointIrrep":
            return func({ unitBuilders, unitFactoryCls, application });
        default:
            throw new Error(`dynamicSubworkflow=${name} not recognized`);
    }
}

function createSubworkflow({
    subworkflowData,
    AppRegistry = ApplicationRegistry,
    modelFactoryCls = ModelFactory,
    methodFactoryCls = MethodFactory,
    subworkflowCls = Subworkflow,
    unitFactoryCls = UnitFactory,
    unitBuilders = builders,
}) {
    const { application, model, method, setSearchText } = createTopLevel({
        subworkflowData,
        AppRegistry,
        modelFactoryCls,
        methodFactoryCls,
    });

    let units = [];
    const { name, units: unitConfigs, config = {}, dynamicSubworkflow = null } = subworkflowData;
    unitConfigs.forEach((_config) => {
        units.push(
            createUnit({
                config: _config,
                application,
                unitBuilders,
                unitFactoryCls,
            }),
        );
    });
    if (dynamicSubworkflow) {
        units = createDynamicUnits({
            dynamicSubworkflow,
            units,
            unitBuilders,
            unitFactoryCls,
            application,
        });
    }

    const { functions = {}, attributes = {}, ...cfg } = config;
    let subworkflow = subworkflowCls.fromArguments(application, model, method, name, units, cfg);
    subworkflow = applyConfig({ obj: subworkflow, config: { functions, attributes } });
    if (setSearchText) subworkflow.model.method.setSearchText(setSearchText);
    return subworkflow;
}

/**
 * @summary Convenience wrapper around createSubworkflow to create by app name and swf name
 * @param appName {String} application name
 * @param swfName {String} subworkflow name (snake_case.yml)
 * @param swArgs {Object} classes for instantiation
 * @returns {*} subworkflow object
 */
function createSubworkflowByName({ appName, swfName, ...swArgs }) {
    const { subworkflows } = allWorkflowData;
    const { [appName]: allSubworkflowData } = subworkflows;
    const { [swfName]: subworkflowData } = allSubworkflowData;
    return createSubworkflow({
        subworkflowData,
        ...swArgs,
    });
}

export { createSubworkflow, createSubworkflowByName };
