import {
    type InMemoryEntityInSet,
    inMemoryEntityInSetMixin,
} from "@mat3ra/code/dist/js/entity/set/InMemoryEntityInSetMixin";
import {
    type OrderedInMemoryEntityInSet,
    orderedEntityInSetMixin,
} from "@mat3ra/code/dist/js/entity/set/ordered/OrderedInMemoryEntityInSetMixin";
import JSONSchemasInterface from "@mat3ra/esse/dist/js/esse/JSONSchemasInterface";
import esseSchemas from "@mat3ra/esse/dist/js/schemas.json";
import type {
    ConditionUnitSchema,
    ErrorUnitSchema,
    ExecutionUnitSchema,
} from "@mat3ra/esse/dist/js/types";
import { Material } from "@mat3ra/made";
import { ApplicationRegistry, WorkflowStandata } from "@mat3ra/standata";
import StandataDriver from "@mat3ra/standata/dist/js/StandataDriver";
import { expect } from "chai";
import type { JSONSchema7 } from "json-schema";
import type { WorkflowRenderContext } from "src/js/Workflow";

import { Subworkflow, UnitFactory, Workflow } from "../../src/js";
import { UnitType } from "../../src/js/enums";
import type { ExecutionUnit } from "../../src/js/units";
import { repairWorkflow } from "../../src/js/utils/repair";
import type { WorkflowEntity } from "../../src/js/Workflow";
import type { WorkflowSchema } from "../../src/js/workflows/types";
import workflowHashes from "../fixtures/workflow_hashes.json";
import { assertNotNull } from "./assertNotNull";

function invalidExecutionUnit(flowchartId: string) {
    return {
        type: UnitType.execution,
        name: "exec",
        flowchartId,
        results: [],
        preProcessors: [],
        postProcessors: [],
        monitors: [],
        input: [],
        context: [],
    };
}

interface OrderedMaterial extends OrderedInMemoryEntityInSet, InMemoryEntityInSet {}

class OrderedMaterial extends Material implements OrderedInMemoryEntityInSet {
    declare static createDefault: () => OrderedMaterial;
}

inMemoryEntityInSetMixin(OrderedMaterial.prototype);
orderedEntityInSetMixin(OrderedMaterial.prototype);

describe("Workflow", () => {
    describe("construction", () => {
        it("creates workflow from default config with name, subworkflows, units, and _id", () => {
            const config = structuredClone(Workflow.defaultConfig);
            const workflow = new Workflow(config);

            expect(workflow.name).to.equal(Workflow.defaultConfig.name);
            expect(workflow.toJSON().subworkflows).to.have.lengthOf(
                Workflow.defaultConfig.subworkflows.length,
            );
            expect(workflow.toJSON().units).to.have.lengthOf(Workflow.defaultConfig.units.length);
            expect(workflow._id).to.be.a("string");
            expect(workflow._id?.length).to.be.above(0);
        });
    });

    describe("ID generation", () => {
        it("assigns different _id to two workflows from default config", () => {
            const w1 = new Workflow(structuredClone(Workflow.defaultConfig));
            const w2 = new Workflow(structuredClone(Workflow.defaultConfig));

            expect(w1._id).to.not.equal(w2._id);
        });
    });

    describe("fromSubworkflow", () => {
        it("creates workflow with one subworkflow and one unit matching subworkflow name", () => {
            const subworkflowConfig = structuredClone(Workflow.defaultConfig.subworkflows[0]);
            const subworkflow = new Subworkflow(subworkflowConfig);
            const workflow = Workflow.fromSubworkflow(subworkflow);

            expect(workflow.toJSON().subworkflows).to.have.lengthOf(1);
            expect(workflow.toJSON().units).to.have.lengthOf(1);
            expect(workflow.name).to.equal(subworkflow.name);
        });
    });

    describe("toJSON", () => {
        it("returns object with name, units, subworkflows, and workflows", () => {
            const config = structuredClone(Workflow.defaultConfig);
            const workflow = new Workflow(config);
            const json = workflow.toJSON();

            expect(json).to.have.property("name", workflow.name);
            expect(json).to.have.property("units").that.is.an("array");
            expect(json).to.have.property("subworkflows").that.is.an("array");
            expect(json).to.have.property("workflows").that.is.an("array");
        });
    });

    describe("getters", () => {
        it("exposes usedApplications, usedApplicationNames, properties, systemName, defaultDescription", () => {
            const config = structuredClone(Workflow.defaultConfig);
            const workflow = new Workflow(config);

            expect(workflow.usedApplications).to.be.an("array");
            expect(workflow.usedApplicationNames).to.be.an("array");
            expect(workflow.getProperties()).to.be.an("array");
            expect(workflow.getSystemName()).to.be.a("string");
            expect(workflow.getSystemName().length).to.be.above(0);
            expect(workflow.getDefaultDescription()).to.be.a("string");
            expect(workflow.getDefaultDescription().length).to.be.above(0);
        });
    });

    describe("addSubworkflow / removeSubworkflow", () => {
        it("adds subworkflows then removes one and updates counts", () => {
            const config = structuredClone(Workflow.defaultConfig);
            const defaultSub = config.subworkflows[0];
            const defaultUnit = config.units[0];
            config.subworkflows = [
                defaultSub,
                { ...defaultSub, _id: "second-sw-id", name: "Second Subworkflow" },
            ];
            config.units = [
                defaultUnit,
                {
                    ...defaultUnit,
                    _id: "second-sw-id",
                    flowchartId: "second-fc-id",
                    name: "Second Subworkflow",
                },
            ];
            const workflow = new Workflow(config);
            const thirdSubworkflow = new Subworkflow({
                ...defaultSub,
                _id: "third-sw-id",
                name: "Third Subworkflow",
            });
            workflow.addSubworkflow(thirdSubworkflow);

            expect(workflow.toJSON().subworkflows).to.have.lengthOf(3);
            expect(workflow.toJSON().units).to.have.lengthOf(3);

            workflow.removeSubworkflow("second-sw-id");

            expect(workflow.toJSON().subworkflows).to.have.lengthOf(2);
            expect(workflow.toJSON().units).to.have.lengthOf(2);
        });
    });

    describe("addUnitType", () => {
        it("adds a subworkflow unit when called with UnitType.subworkflow", () => {
            const config = structuredClone(Workflow.defaultConfig);
            const workflow = new Workflow(config);
            const initialSubworkflows = workflow.toJSON().subworkflows.length;
            const initialUnits = workflow.toJSON().units.length;

            workflow.addUnitType(UnitType.subworkflow);

            expect(workflow.toJSON().subworkflows).to.have.lengthOf(initialSubworkflows + 1);
            expect(workflow.toJSON().units).to.have.lengthOf(initialUnits + 1);
        });
    });

    describe("toggleRelaxation", () => {
        before(() => {
            JSONSchemasInterface.setSchemas(esseSchemas as JSONSchema7[]);
            ApplicationRegistry.setDriver(new StandataDriver());
        });

        it("adds the standata relaxation subworkflow with a schema-complete application", () => {
            const workflow = new Workflow(structuredClone(Workflow.defaultConfig));

            expect(() => workflow.toggleRelaxation()).not.to.throw();

            const relaxation = workflow
                .toJSON()
                .subworkflows.find(
                    (subworkflow) => subworkflow.systemName === "espresso-variable-cell-relaxation",
                );

            if (!relaxation) {
                throw new Error("Expected relaxation subworkflow to be added");
            }

            expect(relaxation.application).to.include({
                build: "GNU",
                name: "espresso",
                shortName: "qe",
                summary: "Quantum ESPRESSO",
                version: "6.3",
            });
        });

        it("inherits the existing subworkflow's application version and build", () => {
            const config = structuredClone(Workflow.defaultConfig);
            Object.assign(config.subworkflows[0].application, { build: "Intel", version: "7.5" });

            const workflow = new Workflow(config);
            workflow.toggleRelaxation();

            const rehydrated = new Workflow(structuredClone(workflow.toJSON()));
            const relaxation = assertNotNull(
                rehydrated.subworkflowInstances.find(
                    (subworkflow) => subworkflow.systemName === "espresso-variable-cell-relaxation",
                ),
            );

            expect(relaxation.application).to.include({ build: "Intel", version: "7.5" });

            const executionUnit = assertNotNull(
                relaxation.unitsInstances.find((unit) => unit.type === UnitType.execution),
            ) as ExecutionUnit;

            expect(executionUnit.application).to.include({ build: "Intel", version: "7.5" });
            expect(executionUnit.executable?.name).to.equal("pw.x");
            expect(executionUnit.flavor?.name).to.equal("pw_vc-relax");
        });

        it("inherits a version and build the registry actually offers", () => {
            const offered = (application: { name: string; version: string; build: string }) =>
                new ApplicationRegistry()
                    .getApplications()
                    .some(
                        (candidate) =>
                            candidate.name === application.name &&
                            candidate.version === application.version &&
                            candidate.build === application.build,
                    );

            const workflow = new Workflow(structuredClone(Workflow.defaultConfig));
            workflow.toggleRelaxation();

            const rehydrated = new Workflow(structuredClone(workflow.toJSON()));
            const relaxation = assertNotNull(
                rehydrated.subworkflowInstances.find(
                    (subworkflow) => subworkflow.systemName === "espresso-variable-cell-relaxation",
                ),
            );

            const executionUnit = assertNotNull(
                relaxation.unitsInstances.find((unit) => unit.type === UnitType.execution),
            ) as ExecutionUnit;

            expect(offered(relaxation.application), "subworkflow application").to.equal(true);
            expect(offered(executionUnit.application), "unit application").to.equal(true);
        });

        it("does not reset the model to the application default when inheriting", () => {
            const config = structuredClone(Workflow.defaultConfig);
            Object.assign(config.subworkflows[0].application, {
                name: "vasp",
                shortName: "vasp",
                summary: "Vienna Ab-initio Simulation Package",
                build: "GNU",
                version: "5.4.4",
            });

            const workflow = new Workflow(config);
            workflow.toggleRelaxation();

            const rehydrated = new Workflow(structuredClone(workflow.toJSON()));
            const relaxation = assertNotNull(
                rehydrated.subworkflowInstances.find(
                    (subworkflow) => subworkflow.systemName === "vasp-variable-cell-relaxation",
                ),
            );

            expect(relaxation.model.method).to.include({
                type: "pseudopotential",
                subtype: "paw",
            });
        });

        it("removes the added relaxation subworkflow on the next toggle", () => {
            const workflow = new Workflow(structuredClone(Workflow.defaultConfig));
            const initialSubworkflowCount = workflow.subworkflowInstances.length;
            const initialUnitCount = workflow.unitInstances.length;

            workflow.toggleRelaxation();

            expect(workflow.hasRelaxation).to.equal(true);
            expect(workflow.subworkflowInstances).to.have.lengthOf(initialSubworkflowCount + 1);
            expect(workflow.unitInstances).to.have.lengthOf(initialUnitCount + 1);

            workflow.toggleRelaxation();

            expect(workflow.hasRelaxation).to.equal(false);
            expect(workflow.subworkflowInstances).to.have.lengthOf(initialSubworkflowCount);
            expect(workflow.unitInstances).to.have.lengthOf(initialUnitCount);
        });
    });

    describe("render", () => {
        before(() => {
            // Context providers resolve JSON Schemas from ESSE at construction time.
            // eslint-disable-next-line global-require, @typescript-eslint/no-var-requires
            JSONSchemasInterface.setSchemas(esseSchemas as JSONSchema7[]);
        });

        it("invokes each subworkflow render with spread context and workflowHasRelaxation", function () {
            this.timeout(10000);
            ApplicationRegistry.setDriver(new StandataDriver());

            const workflows = new WorkflowStandata().getAll().map((standataJson) => {
                return new Workflow(standataJson as unknown as WorkflowSchema);
            });

            expect(
                workflows.length,
                "expected at least one standata workflow from getAll() to construct",
            ).to.be.above(0);

            const material = OrderedMaterial.createDefault();
            const context: WorkflowRenderContext = {
                material,
                materials: [material, material, material],
                jobHasParent: false,
            };

            workflows.forEach((workflow) => {
                workflow.render(context);

                workflow.subworkflowInstances.forEach((subworkflow) => {
                    subworkflow.unitsInstances
                        .filter((unit) => unit.type === UnitType.execution)
                        .forEach((unit) => {
                            expect(unit.renderingContext).to.deep.include({ ...context });
                            expect(unit.renderingContext).to.have.property("methodData");
                            expect(unit.renderingContext).to.have.property("application");
                            expect(unit.renderingContext).to.have.property("subworkflowContext");
                            expect(unit.renderingContext).to.have.property(
                                "workflowHasRelaxation",
                                workflow.hasRelaxation,
                            );
                            expect(unit.context).to.be.deep.equal([]);
                        });
                });
            });
        });
    });

    describe("calculateHash", () => {
        before(() => {
            JSONSchemasInterface.setSchemas(esseSchemas as JSONSchema7[]);
            ApplicationRegistry.setDriver(new StandataDriver());
        });

        // Expected to fail until `toJSON()` / schema clean stops mutating live `_json` used by hashing.
        it("returns the same hash before and after toJSON (standata workflow)", () => {
            const standataWorkflows = new WorkflowStandata().getAll();
            expect(standataWorkflows.length).to.be.above(0);

            const workflow = new Workflow(standataWorkflows[0] as unknown as WorkflowSchema);
            const hashBefore = workflow.calculateHash();
            workflow.toJSON();
            const hashAfter = workflow.calculateHash();

            expect(hashAfter).to.equal(hashBefore);
        });
    });

    describe("Workflow hashing", () => {
        const fixtureFiles = ["band_gap"] as const;
        const bandGapWorkflowName = "Band Gap";

        before(() => {
            JSONSchemasInterface.setSchemas(esseSchemas as JSONSchema7[]);
            ApplicationRegistry.setDriver(new StandataDriver());
        });

        fixtureFiles.forEach((fixtureFile) => {
            it(`calculateHash matches stored hash for ${fixtureFile}`, function () {
                const standata = new WorkflowStandata();
                const workflows = standata.findEntitiesByTags(
                    "espresso",
                    fixtureFile,
                ) as unknown as WorkflowSchema[];
                const workflow = workflows.find((w) => w.name === bandGapWorkflowName);
                if (!workflow) {
                    throw new Error(`Workflow ${bandGapWorkflowName} not found`);
                }
                const wf = new Workflow(workflow);
                const expectedHash = workflowHashes.espresso[fixtureFile].hash;
                if (!expectedHash) {
                    // eslint-disable-next-line no-console
                    console.log(`Hash for espresso/${fixtureFile}: ${wf.calculateHash()}`);
                    this.skip();
                }
                expect(wf.calculateHash()).to.equal(expectedHash);
            });
        });
    });

    describe("repair", () => {
        let defaultApplication: ExecutionUnitSchema["application"];

        before(() => {
            JSONSchemasInterface.setSchemas(esseSchemas as JSONSchema7[]);
        });

        beforeEach(() => {
            ApplicationRegistry.setDriver(new StandataDriver());
            const app = new ApplicationRegistry().getDefaultApplication();
            if (!app) {
                throw new Error(
                    "ApplicationRegistry.getDefaultApplication() returned no application",
                );
            }
            defaultApplication = app;
        });

        it("converts invalid execution unit to error inside subworkflow", () => {
            const workflowConfig = structuredClone(Workflow.defaultConfig);
            const invalid = invalidExecutionUnit("fc-exec");
            workflowConfig.subworkflows[0].units = [invalid as never];

            const result = repairWorkflow(workflowConfig);

            expect(result.subworkflows[0].units[0].type).to.equal(UnitType.error);
            const subworkflowErrorReason = JSON.parse(
                (result.subworkflows[0].units[0] as ErrorUnitSchema).reason,
            );
            expect(subworkflowErrorReason.json).to.deep.equal(invalid);
        });

        it("leaves valid execution unit unchanged inside subworkflow", () => {
            const workflowConfig = structuredClone(Workflow.defaultConfig);
            const unit = UnitFactory.createDefaultSubworkflowUnit(
                "execution",
                defaultApplication,
            ).toJSON();
            workflowConfig.subworkflows[0].units = [unit];

            const result = repairWorkflow(workflowConfig);

            expect(result.subworkflows[0].units[0]).to.deep.equal(unit);
        });

        it("repairs legacy condition unit missing then and else via repairWorkflow", () => {
            const workflowConfig = structuredClone(Workflow.defaultConfig);
            const legacyCondition = {
                name: "condition",
                type: UnitType.condition,
                input: [],
                results: [],
                preProcessors: [],
                postProcessors: [],
                statement: "true",
                maxOccurrences: 100,
                application: defaultApplication,
                flowchartId: "e71f01a98db152d645e787b8",
                monitors: [],
                head: false,
                next: "03b383bf35936d1d54d015ad",
            };

            workflowConfig.subworkflows[0].units = [legacyCondition as never];

            const result = repairWorkflow(workflowConfig);
            const condition = result.subworkflows[0].units[0] as ConditionUnitSchema;

            expect(condition.type).to.equal(UnitType.condition);
            expect(condition.then).to.equal("");
            expect(condition.else).to.equal("");
            expect(() => new Workflow(result).validate()).to.not.throw();
        });

        it("repairs only invalid execution units inside subworkflow", () => {
            const workflowConfig = structuredClone(Workflow.defaultConfig);
            const validUnit = UnitFactory.createDefaultSubworkflowUnit("assignment").toJSON();
            const invalidExecution = invalidExecutionUnit("fc-bad");
            workflowConfig.subworkflows[0].units = [validUnit, invalidExecution as never];

            const result = repairWorkflow(workflowConfig);

            expect(result.subworkflows[0].units[0]).to.deep.equal(validUnit);
            expect(result.subworkflows[0].units[1].type).to.equal(UnitType.error);
        });

        it("repairs invalid execution units in subworkflows and nested workflows", () => {
            const workflowConfig = structuredClone(Workflow.defaultConfig);
            const invalidExecution = invalidExecutionUnit("fc-nested");

            workflowConfig.subworkflows[0].units = [invalidExecution as never];

            const nestedSubworkflow = structuredClone(workflowConfig.subworkflows[0]);
            nestedSubworkflow._id = "nested-subworkflow-id";
            nestedSubworkflow.units = [invalidExecution as never];

            workflowConfig.workflows = [
                {
                    ...structuredClone(workflowConfig),
                    _id: "nested-workflow-id",
                    subworkflows: [nestedSubworkflow],
                    workflows: [],
                },
            ];

            const result = repairWorkflow(workflowConfig);

            expect(result.subworkflows[0].units[0].type).to.equal(UnitType.error);
            expect(result.workflows[0].subworkflows[0].units[0].type).to.equal(UnitType.error);
        });

        it("allows hydration after repair for converted execution units", () => {
            const workflowConfig = structuredClone(Workflow.defaultConfig);
            workflowConfig.subworkflows[0].units = [invalidExecutionUnit("fc-hydrate") as never];

            const document = repairWorkflow(workflowConfig);

            expect(() => new Workflow(document)).to.not.throw();
        });

        it("repairs legacy subworkflow missing model without throwing", () => {
            const legacyShellWorkflow = {
                name: "Bash test",
                subworkflows: [
                    {
                        name: "bash subworkflow",
                        application: {
                            name: "shell",
                            version: "4.2.46",
                            build: "Default",
                            isDefault: true,
                            summary: "Shell Script",
                            shortName: "sh",
                        },
                        units: [
                            {
                                type: UnitType.execution,
                                name: "test workflow",
                                head: true,
                                application: {
                                    name: "shell",
                                    version: "4.2.46",
                                    build: "Default",
                                    isDefault: true,
                                    shortName: "sh",
                                    summary: "Shell Script",
                                },
                                results: [],
                                monitors: [{ name: "standard_output" }],
                                preProcessors: [],
                                postProcessors: [],
                                input: [],
                                context: [],
                            },
                        ],
                        isDraft: false,
                        properties: [],
                    },
                ],
                units: [
                    {
                        name: "bash subworkflow",
                        type: UnitType.subworkflow,
                        head: true,
                        monitors: [],
                        results: [],
                        preProcessors: [],
                        postProcessors: [],
                    },
                ],
                workflows: [],
                properties: [],
            } as unknown as WorkflowSchema;

            const result = repairWorkflow(legacyShellWorkflow);

            expect(result.subworkflows).to.have.lengthOf(0);
            expect(result.units[0].type).to.equal(UnitType.error);
            const legacyErrorReason = JSON.parse((result.units[0] as ErrorUnitSchema).reason);
            expect(legacyErrorReason.error).to.be.an("array").that.is.not.empty;
            expect(() => new Workflow(result)).to.not.throw();
        });

        it("does not store JSON schema in error unit reason", () => {
            const workflowConfig = structuredClone(Workflow.defaultConfig);
            const invalid = invalidExecutionUnit("fc-1");
            workflowConfig.subworkflows[0].units = [invalid as never];

            const result = repairWorkflow(workflowConfig);
            const errorUnit = result.subworkflows[0].units[0] as ErrorUnitSchema;
            const reason = JSON.parse(errorUnit.reason);

            expect(reason).to.have.property("error");
            expect(reason).to.have.property("json");
            expect(reason).to.not.have.property("schema");
        });

        it("converts invalid subworkflow unit to error and drops subworkflow", () => {
            const workflowConfig = structuredClone(Workflow.defaultConfig);
            const originalUnit = workflowConfig.units[0];
            const subworkflowId = workflowConfig.subworkflows[0]._id;

            workflowConfig.subworkflows[0] = {
                ...workflowConfig.subworkflows[0],
                model: {
                    ...workflowConfig.subworkflows[0].model,
                    type: "not-a-valid-model-type" as never,
                },
            };

            const result = repairWorkflow(workflowConfig);

            expect(result.subworkflows).to.have.lengthOf(0);
            expect(result.units).to.have.lengthOf(1);
            expect(result.units[0].type).to.equal(UnitType.error);
            expect(result.units[0]._id).to.equal(subworkflowId);

            const errorUnit = result.units[0] as ErrorUnitSchema;
            const validationReason = JSON.parse(errorUnit.reason);
            expect(validationReason.error).to.be.an("array").that.is.not.empty;
            expect(validationReason.error.some((e: { instancePath?: string }) => e.instancePath)).to
                .be.true;
            expect(validationReason.json).to.deep.equal({
                unit: originalUnit,
                subworkflow: workflowConfig.subworkflows[0],
            });
            expect(() => new Workflow(result)).to.not.throw();
        });
    });

    describe("generic schema wrapper", () => {
        type WiderWorkflowSchema = WorkflowEntity & { webappOnly?: string };

        class WiderWorkflow extends Workflow<WiderWorkflowSchema> {}

        it("allows subclasses to widen _json typing and storage", () => {
            const workflow = new WiderWorkflow(structuredClone(Workflow.defaultConfig));

            workflow._json.webappOnly = "webapp-value";
            expect(workflow._json.webappOnly).to.equal("webapp-value");

            expect(workflow.toJSON().name).to.equal(Workflow.defaultConfig.name);
            expect(workflow._json.webappOnly).to.equal("webapp-value");

            workflow._json.webappOnly = "updated";
            expect(workflow._json.webappOnly).to.equal("updated");
        });
    });
});
