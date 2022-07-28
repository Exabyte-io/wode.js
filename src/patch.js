import { Template } from "@exabyte-io/ade.js";

import { ContextProviderRegistry } from "./context/registry";

// We patch the static providerRegistry here so that
// Template has all context providers available
// to it when creating workflows. It is then re-exported
// from WoDe for use downstream.
Template.providerRegistry = ContextProviderRegistry;

export { Template };
