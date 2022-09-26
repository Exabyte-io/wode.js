import { UNIT_TYPES } from "../enums";

export const ConvergenceMixin = (superclass) =>
    class extends superclass {
        addConvergence({
            parameter,
            parameterInitial,
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

            // Replace kgrid to be ready for convergence
            // TODO: kgrid should be abstracted and selected by user
            unitForConvergence.updateContext({
                kgrid: {
                    dimensions: [`{{${parameter}}}`, `{{${parameter}}}`, `{{${parameter}}}`],
                    shifts: [0, 0, 0],
                },
                isKgridEdited: true,
                usesJinjaVariable: true,
            });

            const prevResult = "prev_result";

            // Assignment with result's initial value
            const prevResultInit = this._UnitFactory.create({
                type: UNIT_TYPES.assignment,
                head: true,
                operand: prevResult,
                value: resultInitial,
            });

            // Assignment with initial value of convergence parameter
            const paramInit = this._UnitFactory.create({
                type: UNIT_TYPES.assignment,
                operand: parameter,
                value: parameterInitial,
            });

            // Assignment for storing iteration result: extracts 'result' from convergence unit scope
            const storePrevResult = this._UnitFactory.create({
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
                type: UNIT_TYPES.assignment,
                input: [],
                operand: parameter,
                value: `${parameter} + 1`,
                next: unitForConvergence.flowchartId,
            });

            // Final step of convergence
            const exit = this._UnitFactory.create({
                type: UNIT_TYPES.assignment,
                name: "exit",
                input: [],
                operand: parameter,
                value: `${parameter} + 0`,
            });

            // Final step of convergence
            const storeResult = this._UnitFactory.create({
                type: UNIT_TYPES.assignment,
                input: [
                    {
                        scope: unitForConvergence.flowchartId,
                        name: result,
                    },
                ],
                operand: result,
                value: result,
            });

            // Convergence condition unit
            const conditionUnit = this._UnitFactory.create({
                type: UNIT_TYPES.condition,
                statement: `${condition} ${operator} ${tolerance}`,
                then: exit.flowchartId,
                else: storePrevResult.flowchartId,
                maxOccurrences,
                next: storePrevResult.flowchartId,
            });

            this.addUnit(prevResultInit, true);
            this.addUnit(paramInit, true);
            this.addUnit(storeResult);
            this.addUnit(conditionUnit);
            this.addUnit(storePrevResult);
            this.addUnit(nextStep);
            this.addUnit(exit);

            // `addUnit` adjusts the `next` field, hence the below.
            nextStep.next = unitForConvergence.flowchartId;
        }
    };
