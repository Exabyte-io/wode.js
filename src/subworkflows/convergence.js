import merge from "lodash/merge";

import { UNIT_TAGS, UNIT_TYPES } from "../enums";
import { createConvergenceParameter } from "./convergence/factory";

export const ConvergenceMixin = (superclass) =>
    class extends superclass {
        get convergenceParam() {
            return this.findUnitWithTag("hasConvergenceParam")?.operand || undefined;
        }

        get convergenceResult() {
            return this.findUnitWithTag("hasConvergenceResult")?.operand || undefined;
        }

        convergenceSeries(scopeTrack) {
            if (!this.hasConvergence || !scopeTrack?.length) return [];
            let lastResult;
            const series = scopeTrack
                .map((scopeItem, i) => ({
                    x: i,
                    param: scopeItem.scope?.global[this.convergenceParam],
                    y: scopeItem.scope?.global[this.convergenceResult],
                }))
                .filter(({ y }) => {
                    const isNewResult = y !== undefined && y !== lastResult;
                    lastResult = y;
                    return isNewResult;
                });
            return series.map((item, i) => {
                return {
                    x: i + 1,
                    param: item.param,
                    y: item.y,
                };
            });
        }

        addConvergence({
            parameter,
            parameterInitial,
            parameterIncrement,
            result,
            resultInitial,
            condition,
            operator,
            tolerance,
            maxOccurrences,
        }) {
            // RF: added TODO comments for future improvements

            const { units } = this;
            // Find unit to converge: should contain passed result in its results list
            // TODO: make user to select unit for convergence explicitly
            const unitForConvergence = units.find((x) =>
                x.resultNames.find((name) => name === result),
            );

            if (!unitForConvergence) {
                // eslint-disable-next-line no-undef
                sAlert.error(
                    `Subworkflow does not contain unit with '${result}' as extracted property.`,
                );
                throw new Error("There is no result to converge");
            }

            // initialize parameter
            const param = createConvergenceParameter({
                name: parameter,
                initialValue: parameterInitial,
                increment: parameterIncrement,
            });

            // Replace kgrid to be ready for convergence
            // TODO: kgrid should be abstracted and selected by user
            const providers = unitForConvergence.importantSettingsProviders;
            const gridProvider = providers.find((p) => p.name === "kgrid" || p.name === "qgrid");
            let mergedContext = param.unitContext;
            if (gridProvider) {
                mergedContext = merge(gridProvider.yieldData(), param.unitContext);
                gridProvider.setData(mergedContext);
                gridProvider.setIsEdited(true);
            }
            unitForConvergence.updateContext(mergedContext);

            const prevResult = "prev_result";
            const iteration = "iteration";

            // Assignment with result's initial value
            const prevResultInit = this._UnitFactory.create({
                name: "init result",
                type: UNIT_TYPES.assignment,
                head: true,
                operand: prevResult,
                value: resultInitial,
            });

            // Assignment with initial value of convergence parameter
            const paramInit = this._UnitFactory.create({
                name: "init parameter",
                type: UNIT_TYPES.assignment,
                operand: param.name,
                value: param.initialValue,
                tags: [UNIT_TAGS.hasConvergenceParam],
            });

            // Assignment with initial value of iteration counter
            const iterInit = this._UnitFactory.create({
                name: "init counter",
                type: UNIT_TYPES.assignment,
                operand: iteration,
                value: 1,
            });

            // Assignment for storing iteration result: extracts 'result' from convergence unit scope
            const storePrevResult = this._UnitFactory.create({
                name: "store result",
                type: UNIT_TYPES.assignment,
                input: [
                    {
                        scope: unitForConvergence.flowchartId,
                        name: result,
                    },
                ],
                operand: prevResult,
                value: result,
            });

            // Assignment for convergence param increase
            const nextStep = this._UnitFactory.create({
                name: "update parameter",
                type: UNIT_TYPES.assignment,
                input: param.useVariablesFromUnitContext(unitForConvergence.flowchartId),
                operand: param.name,
                value: param.increment,
                next: unitForConvergence.flowchartId,
            });

            // Final step of convergence
            const exit = this._UnitFactory.create({
                type: UNIT_TYPES.assignment,
                name: "exit",
                input: [],
                operand: param.name,
                value: param.finalValue,
            });

            // Final step of convergence
            const storeResult = this._UnitFactory.create({
                name: "update result",
                type: UNIT_TYPES.assignment,
                input: [
                    {
                        scope: unitForConvergence.flowchartId,
                        name: result,
                    },
                ],
                operand: result,
                value: result,
                tags: [UNIT_TAGS.hasConvergenceResult],
            });

            // Assign next iteration value
            const nextIter = this._UnitFactory.create({
                name: "update counter",
                type: UNIT_TYPES.assignment,
                input: [],
                operand: iteration,
                value: `${iteration} + 1`,
            });

            // Convergence condition unit
            const conditionUnit = this._UnitFactory.create({
                name: "check convergence",
                type: UNIT_TYPES.condition,
                statement: `${condition} ${operator} ${tolerance}`,
                then: exit.flowchartId,
                else: storePrevResult.flowchartId,
                maxOccurrences,
                next: storePrevResult.flowchartId,
            });

            this.addUnit(paramInit, 0);
            this.addUnit(prevResultInit, 1);
            this.addUnit(iterInit, 2);
            this.addUnit(storeResult);
            this.addUnit(conditionUnit);
            this.addUnit(storePrevResult);
            this.addUnit(nextIter);
            this.addUnit(nextStep);
            this.addUnit(exit);

            // `addUnit` adjusts the `next` field, hence the below.
            nextStep.next = unitForConvergence.flowchartId;
        }
    };
