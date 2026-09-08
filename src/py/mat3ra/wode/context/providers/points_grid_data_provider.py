import math
from typing import Any, Dict, List, Optional

from mat3ra.esse.models.context_providers_directory.points_grid_data_provider import (
    GridMetricType,
    PointsGridDataProviderSchema,
)
from mat3ra.made.material import Material
from pydantic import Field, model_validator

from .base.context_provider import ContextProvider

DEFAULT_KPPRA = -1

SCHEMA_FIELD_NAMES = set(PointsGridDataProviderSchema.model_fields)


# TODO: GlobalSetting for default KPPRA value
class PointsGridDataProvider(PointsGridDataProviderSchema, ContextProvider):
    """Context provider for k-point/q-point grid configuration."""

    name: str = Field(default="kgrid")
    divisor: int = Field(default=1)
    dimensions: List[int] = Field(default_factory=lambda: [1, 1, 1], min_length=3, max_length=3)
    shifts: List[float] = Field(default_factory=lambda: [0.0, 0.0, 0.0])
    gridMetricType: GridMetricType = Field(default=GridMetricType.KPPRA)
    gridMetricValue: float = Field(default=DEFAULT_KPPRA)
    preferGridMetric: bool = Field(default=False)
    material: Optional[Material] = Field(default=None, exclude=True)

    @model_validator(mode="after")
    def _derive_grid_metric_and_ratios(self) -> "PointsGridDataProvider":
        if "dimensions" in self.model_fields_set and "gridMetricValue" not in self.model_fields_set:
            self.gridMetricValue = self.calculate_grid_metric(self.gridMetricType, self.dimensions)
        if self.reciprocalVectorRatios is None and self.material is not None:
            # rounded like JS, so job context matches what the UI writes
            ratios = self.material.lattice.reciprocal_vector_ratios
            self.reciprocalVectorRatios = [round(float(r), 3) for r in ratios]
        return self

    _MATERIAL_REQUIRED = "KPPRA and reciprocal vector ratios come from material=<mat3ra.made.Material>"

    def get_number_of_atoms(self) -> int:
        if self.material is None:
            raise ValueError(self._MATERIAL_REQUIRED)
        return self.material.basis.number_of_atoms

    @property
    def is_edited_key(self) -> str:
        return "isKgridEdited"

    @property
    def default_data(self) -> Dict[str, Any]:
        return self.model_dump(by_alias=True, exclude_none=True, include=SCHEMA_FIELD_NAMES)

    def get_reciprocal_vector_ratios(self, context: Optional[Dict[str, Any]] = None) -> Optional[List[float]]:
        effective_data = self._get_effective_data(context)
        if isinstance(effective_data, dict) and "reciprocalVectorRatios" in effective_data:
            return effective_data["reciprocalVectorRatios"]
        return self.reciprocal_vector_ratios

    def build_data(
        self,
        *,
        dimensions: Optional[List[Any]] = None,
        shifts: Optional[List[float]] = None,
        reciprocal_vector_ratios: Optional[List[float]] = None,
        grid_metric_type: Optional[str] = None,
        grid_metric_value: Optional[float] = None,
        prefer_grid_metric: Optional[bool] = None,
        divisor: Optional[int] = None,
    ) -> Dict[str, Any]:
        data = dict(self.default_data)
        overrides = {
            "dimensions": dimensions,
            "shifts": shifts,
            "reciprocalVectorRatios": reciprocal_vector_ratios,
            "gridMetricType": grid_metric_type,
            "gridMetricValue": grid_metric_value,
            "preferGridMetric": prefer_grid_metric,
            "divisor": divisor,
        }
        data.update({key: value for key, value in overrides.items() if value is not None})
        return data

    def yield_data_with_overrides(
        self,
        *,
        dimensions: Optional[List[Any]] = None,
        shifts: Optional[List[float]] = None,
        reciprocal_vector_ratios: Optional[List[float]] = None,
        grid_metric_type: Optional[str] = None,
        grid_metric_value: Optional[float] = None,
        prefer_grid_metric: Optional[bool] = None,
        divisor: Optional[int] = None,
        is_using_jinja_variables: bool = False,
    ) -> Dict[str, Any]:
        context = self.yield_data(
            context={
                self.name_str: self.build_data(
                    dimensions=dimensions,
                    shifts=shifts,
                    reciprocal_vector_ratios=reciprocal_vector_ratios,
                    grid_metric_type=grid_metric_type,
                    grid_metric_value=grid_metric_value,
                    prefer_grid_metric=prefer_grid_metric,
                    divisor=divisor,
                ),
                self.is_edited_key: True,
            }
        )
        if is_using_jinja_variables:
            context["isUsingJinjaVariables"] = True
        return context

    # TODO: add a test to verify context and templates are the same as from JS implementation
    def get_default_grid_metric_value(self, metric: str) -> float:
        raise NotImplementedError

    def calculate_dimensions(
        self, grid_metric_type: GridMetricType, grid_metric_value: float, units: str = "angstrom"
    ) -> List[int]:
        raise NotImplementedError

    def calculate_grid_metric(self, grid_metric_type: GridMetricType, dimensions: List[int]) -> float:
        if grid_metric_type == GridMetricType.KPPRA:
            return math.prod(dimensions) * self.get_number_of_atoms()
        raise NotImplementedError(f"calculate_grid_metric not implemented for {grid_metric_type}")

    def transform_data(self, data: Dict[str, Any]) -> Dict[str, Any]:
        raise NotImplementedError
