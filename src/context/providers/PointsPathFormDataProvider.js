/* eslint-disable max-classes-per-file */
import { Application } from "@exabyte-io/ade.js";
import {
    ApplicationContextMixin,
    JSONSchemaFormDataProvider,
    MaterialContextMixin,
} from "@exabyte-io/code.js/dist/context";
import { math as codeJSMath } from "@exabyte-io/code.js/dist/math";
import { Made } from "@exabyte-io/made.js";
import { mix } from "mixwith";
import s from "underscore.string";

const defaultPoint = "Г";
const defaultSteps = 10;

export class PointsPathFormDataProvider extends mix(JSONSchemaFormDataProvider).with(
    ApplicationContextMixin,
    MaterialContextMixin,
) {
    static Material = Made.Material;

    static Application = Application;

    constructor(config) {
        super(config);
        this.reciprocalLattice = new Made.ReciprocalLattice(this.material.lattice);
        this.symmetryPoints = this.symmetryPointsFromMaterial;
    }

    get isEditedIsSetToFalseOnMaterialUpdate() {
        return this.isMaterialUpdated || this.isMaterialCreatedDefault;
    }

    get defaultData() {
        return this.reciprocalLattice.defaultKpointPath;
    }

    get symmetryPointsFromMaterial() {
        return this.reciprocalLattice.symmetryPoints;
    }

    get jsonSchema() {
        // no need to pass context to get symmetry points on client
        const points = [].concat(this.symmetryPoints).map((x) => x.point);
        return {
            $schema: "http://json-schema.org/draft-04/schema#",
            title: " ",
            description: "path in reciprocal space",
            type: "array",
            items: {
                type: "object",
                properties: {
                    point: {
                        type: "string",
                        default: defaultPoint,
                        enum: points,
                    },
                    steps: {
                        type: "integer",
                        default: defaultSteps,
                    },
                },
            },
            minItems: 1,
        };
    }

    // eslint-disable-next-line class-methods-use-this
    get uiSchema() {
        return {
            items: {},
        };
    }

    get uiSchemaStyled() {
        return {
            items: {
                point: this.defaultFieldStyles,
                steps: this.defaultFieldStyles,
            },
        };
    }

    get fields() {
        const hasRequiredFn = typeof this.material.getBrillouinZoneImageComponent === "function";
        if (!hasRequiredFn) {
            console.log(
                "PointsPathFormDataProvider: Material class has no function" +
                    " 'getBrillouinZoneImageComponent'! Returning empty Object instead.",
            );
            return {};
        }
        return {
            // eslint-disable-next-line no-unused-vars
            TitleField: ({ title, required }) =>
                this.material.getBrillouinZoneImageComponent(title),
        };
    }

    get useExplicitPath() {
        return this.application.name === "vasp";
    }

    // override yieldData to avoid storing explicit path in saved context
    yieldDataForRendering() {
        return this.yieldData(this.useExplicitPath);
    }

    transformData(path = [], useExplicitPath = false) {
        const rawData = path.map((p) => {
            const point = this.symmetryPoints.find((sp) => sp.point === p.point);
            return { ...p, coordinates: point.coordinates };
        });
        const processedData = useExplicitPath ? this._convertToExplicitPath(rawData) : rawData;
        // make coordinates into string and add formatting
        return processedData.map((p) => {
            const coordinates = this.is2PIBA
                ? this.get2PIBACoordinates(p.coordinates)
                : p.coordinates;
            p.coordinates = coordinates.map((c) => s.sprintf("%14.9f", c));
            return p;
        });
    }

    get2PIBACoordinates(point) {
        return this.reciprocalLattice.getCartesianCoordinates(point);
    }

    // Initially, path contains symmetry points with steps counts.
    // This function explicitly calculates each point between symmetry points by step counts.
    // eslint-disable-next-line class-methods-use-this
    _convertToExplicitPath(path) {
        const points = [];
        for (let i = 0; i < path.length - 1; i++) {
            const startPoint = path[i];
            const endPoint = path[i + 1];
            const middlePoints = codeJSMath.calculateSegmentsBetweenPoints3D(
                startPoint.coordinates,
                endPoint.coordinates,
                startPoint.steps,
            );
            points.push(startPoint.coordinates);
            points.push(...middlePoints);
            // Include endPoint into path for the last section, otherwise it will be included by next loop iteration
            if (path.length - 2 === i) points.push(endPoint.coordinates);
        }
        return points.map((x) => {
            return {
                coordinates: x,
                steps: 1,
            };
        });
    }
}

export class ExplicitPointsPathFormDataProvider extends PointsPathFormDataProvider {
    // eslint-disable-next-line class-methods-use-this
    get useExplicitPath() {
        return true;
    }
}

export class ExplicitPointsPath2PIBAFormDataProvider extends ExplicitPointsPathFormDataProvider {
    // eslint-disable-next-line class-methods-use-this
    get is2PIBA() {
        return true;
    }
}
